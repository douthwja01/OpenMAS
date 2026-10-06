# Agents and objects

Every participant in an OpenMAS scene is a MATLAB handle class under `objects/`. The simulator does not call into private decision code. It calls a small public surface: construct the object, `setup` once, then `main` once per step. The scene moves by integrating `GLOBAL.velocity`:

```text
position = position + velocity * dt
```

`main` therefore publishes a new velocity and quaternion before it returns. The position stored on the object is the agent's own copy of that integral. The position used for sensing, collision, and `DATA.globalTrajectories` is the copy kept on `META.OBJECTS`, advanced by `OMAS_updateGlobalStates` at the start of the next step.

See [Simulation guide](README.md) for the loop that calls these methods, the observation packet, and `OMAS_initialise`.

## Class hierarchy

```text
objectDefinition          handle base: identity, geometry, GLOBAL pose, 6-DoF default motion
├── agent                 sensing, memory, waypoints, kinematic limits
│   └── agent_2D          planar state, memory, and control
│       └── agent_2D_*    VO, RVO, HRVO, ORCA, interval, formation, ...
├── obstacle              passive body, spherical hit box
│   ├── obstacle_cuboid
│   └── obstacle_spheroid
└── waypoint              passive goal with an ownership list
```

`agent` also inherits `agent_tools`, which holds memory, waypoint selection, and the global-update helpers shared by every active object. Concrete algorithms (`agent_VO`, `agent_2D_ORCA`, `quadcopter`, `ARdrone`, `fixedWing`, ...) override `main` and, when the motion model is different, the dynamics or `GlobalUpdate` method.

Put new agent classes in `objects/agents/` and add the objects tree with `addpath(genpath('objects'))`. The constructor walks `superclasses` when it looks for an STL mesh under `objects/models/` named after the class or one of its parents.

## What the simulator reads and writes

`GLOBAL` is a private struct on `objectDefinition`. Use `GetGLOBAL` and `SetGLOBAL` rather than touching it directly.

```matlab
pose = agent.GetGLOBAL();                  % whole struct
p    = agent.GetGLOBAL('position');        % one field
agent.SetGLOBAL('position', [0; 0; 0]);    % one field
agent.SetGLOBAL(pose);                     % replace the struct
```

Fields the simulator requires before the run, and reads again after every `main`:

| Field | Shape | Role |
| --- | --- | --- |
| `type` | `OMAS_objectType` | Selects the update branch: agent, obstacle, waypoint, or misc. |
| `hitBoxType` | `OMAS_hitBoxType` | Collision volume. `none` skips warning and collision tests. |
| `position` | `3x1` | Object-side global position, East-North-Up, metres. The recorded trajectory integrates velocity separately. |
| `velocity` | `3x1` | Global velocity, metres per second. This is the value the simulator integrates. |
| `quaternion` | `4x1` | Attitude of the body axes in the global frame. Identity is `[1;0;0;0]`. Rebuilt into `R` on the next step. |
| `radius` | scalar | Characteristic size. Also scales STL geometry. |
| `detectionRadius` | scalar | Sensing horizon in metres. `inf` sees the whole scene. Agents only. |
| `idleStatus` | logical | `true` when the agent has finished its task. All-idle ends the run after `idleTimeOut`. |
| `is3D` | logical | `true` keeps 3-vectors in observations. `agent_2D` sets this `false`. |
| `colour` | `1x3` | RGB used in figures. |
| `symbol` | char / string | Marker used in plan plots. |

`radius` and `detectionRadius` have setters on the base classes. Assigning `this.radius = 1` updates `GLOBAL.radius` and rescales any loaded mesh. Assigning `this.detectionRadius = 30` updates `GLOBAL.detectionRadius` and `SENSORS.range`.

`objectID` is allocated automatically and is stable for the life of the MATLAB session. `name` defaults to a Greek-letter tag (`alpha001`, ...) when you do not pass `'name'`. Names cannot contain `\ / : * ? " < > |` because they are used in file names.

## Lifecycle

### 1. Construction

Name-value pairs override public properties. The base constructor calls `ApplyUserOverrides`, which matches those names against properties.

```matlab
a = agent_example('radius', 0.5, 'detectionRadius', 40, 'name', 'alpha', 'v_max', 3);
```

A subclass constructor should call the superclass first, set its own defaults, then apply user overrides again so the caller's pairs win:

```matlab
function [this] = myAgent(varargin)
    this@agent(varargin);
    this.radius = 0.5;
    this.detectionRadius = 50;
    this.v_nominal = 2;
    this.v_max = 4;
    this = this.ApplyUserOverrides(varargin);
end
```

At this moment the object is still at the origin with zero velocity and an identity quaternion. Scenario code places it afterwards.

### 2. Scenario placement

`GetScenario_*` or your own script writes the initial world state:

```matlab
agent.SetGLOBAL('position',   [x; y; z]);
agent.SetGLOBAL('velocity',   [vx; vy; vz]);
agent.SetGLOBAL('quaternion', [qw; qx; qy; qz]);
```

`OMAS_initialise` checks those three vectors, builds the rotation matrix with `quat2rotm`, and converts velocity and ZYX Euler angles into the body frame.

### 3. `setup` (once)

```matlab
this = this.setup(localXYZvelocity, localXYZrotations);
```

`agent.setup` dispatches on `Is3D()`:

- 3D state is `12x1`: `[x; y; z; phi; theta; psi; dx; dy; dz; dphi; dtheta; dpsi]`.
- 2D state is `6x1`: `[x; y; psi; dx; dy; dpsi]`.

`setup_3DVelocities` stores the body-frame initial velocity in `localState(7:9)` and the roll/pitch part of the Euler vector in `localState(4:5)`. Positions inside `localState` stay zero; the world position lives on `GLOBAL`. Override `setup` when the state vector has a different layout, and still call `SetGLOBAL('priorState', this.localState)` so later rate estimates have a previous sample.

### 4. `main` (every step)

This is the only method the time loop calls on an agent.

```matlab
function [this] = main(this, ENV, varargin)
    % ENV.dt, ENV.currentTime, ENV.currentStep, ENV.timeVector, ...
    observations = varargin{1};   % struct array, or []

    [this, obstacleSet, agentSet] = this.GetAgentUpdate(ENV, observations);

    heading = this.GetTargetHeading();
    desiredVelocity = heading * this.v_nominal;   % body-frame 3-vector
    this = this.Controller(ENV.dt, desiredVelocity);
end
```

`GetAgentUpdate` is the standard way to absorb a packet:

1. Return immediately when the packet is empty.
2. Run each detection through `SensorModel`.
3. Append the measurement to `MEMORY` with `UpdateMemoryFromObject`.
4. Sort memory by `priority`.
5. Split memory into `agentSet`, `waypointSet`, and `obstacleSet` using `OMAS_objectType`.
6. Call `UpdateTargetWaypoint`, which tracks the highest-priority associated waypoint and sets `idleStatus` when none remain.

`SensorModel` is the override point for imperfect sensing. The default copies camera angles (`heading`, `elevation`, `width`) and range, then adds Gaussian noise from `SENSORS`. The default sigmas are zero, so the packet is passed through unchanged. `GetCustomSensorParameters` is a noisier preset (`sigma_position = 0.5`, `sigma_velocity = 0.1`, and a one-pixel camera sigma).

The base `agent.main` does the same update, aims at the current waypoint at `v_nominal`, and calls `Controller`. Replace `main` in a subclass and keep the two obligations: consume `ENV` and the packet, and publish a new global pose before returning.

## Publishing motion

`localState` is the agent's private estimate. `OMAS_updateGlobalStates` runs before `main` and reads only velocity, quaternion, and idle status. The velocity published by this call of `main` is applied as `position + velocity * dt` at the beginning of the next step, and the quaternion becomes the body frame for the next observation packet.

These methods integrate a local decision into the object's `GLOBAL` pose so that the next read sees a consistent velocity. All of them finish by calling `GlobalUpdate_direct(position, velocity, quaternion)`, which also stores `R` and the previous `localState` as `priorState`. The object's position is updated with the previous global velocity times `dt`, matching the simulator's integral when both sides start from the same state.

| Method | State it expects | What it integrates |
| --- | --- | --- |
| `GlobalUpdate(dt, eulerState)` | `[x y z phi theta psi]` or `[x y psi]` | Finite-difference linear and angular rates from `priorState`, then rotates them into the global frame. |
| `GlobalUpdate_3DVelocities(dt, eulerState)` | 12-state vector with rates in the last six entries | Uses `localState(7:12)` as body rates directly. This is what `Controller` calls. |
| `GlobalUpdate_fixedFrame(dt, eulerState)` | Position states | Updates translation and keeps the current quaternion. |
| `GlobalUpdate_direct(p, v, q)` | Global `3x1`, `3x1`, `4x1` | Writes the world state with no integration. Use this for a custom dynamic model that already produces global quantities. |

`Controller(dt, desiredVelocity)` is the built-in kinematic tracker. `desiredVelocity` is a body-frame 3-vector. The controller:

1. Converts the direction into yaw and pitch rates relative to body `x`.
2. Clamps speed and turn rate with `ApplyKinematicContraints` (`v_max`, `w_max`, and the `DYNAMICS` acceleration limits).
3. Holds position when `IsIdle()` is true.
4. Steps `localState` with `SimpleDynamics` and publishes it through `GlobalUpdate_3DVelocities`.

`Controller_PID` is the same idea with proportional feedback on speed error and heading error. Default limits are `v_nominal = 2`, `v_max = 4`, `w_max = 0.5` rad/s, and infinite linear and angular acceleration. Set `v_max` or `w_max` and the matching `DYNAMICS` field updates with them.

`SimpleDynamics(state, linearVelocity, angularVelocity)` returns the state derivative for a single integrator: the linear velocity is written into the position-rate slots and the angular velocity into the attitude-rate slots. It does not rotate anything. Rotation into East-North-Up happens inside `GlobalUpdate_*`.

A custom dynamic model (quadcopter, fixed wing, ARDrone) should integrate its own state inside `main`, then call the global update that matches that state. `quadcopter` and `ARdrone` override `GlobalUpdate` for this reason.

## Memory

`agent.MEMORY` is a struct array, one entry per observed object, created by `CreateMEMORY`. Cartesian fields are 3-wide for a 3D agent and 2-wide for `agent_2D`. Each measurement is a `circularBuffer` of length `maxSamples` (default 1, so only the latest sample is kept).

| Field | Meaning |
| --- | --- |
| `name`, `objectID`, `type` | Identity of the observed object. |
| `sampleNum`, `time` | How many hits have been stored, and when. |
| `position`, `velocity`, `radius` | Body-frame Cartesian measurement. |
| `range`, `heading`, `elevation`, `width` | Spherical measurement after `SensorModel`. |
| `geometry` | Relative mesh, when one was sent. |
| `priority` | Waypoint priority for this agent. Empty for other types. |

`GetLastMeasurementByID(id, field)` reads the newest sample. `UpdateMemoryOrderByField('priority')` sorts the array. Entries are not removed when an object leaves view; a later detection with the same `objectID` appends another sample.

After `GetAgentUpdate`, the three sets are views of that memory, selected by `type`. Their `.position` vectors are in the agent body frame, so a relative-velocity obstacle algorithm can use them directly.

## Waypoints and idle

A `waypoint` is an `objectDefinition` with type `waypoint`, a spherical hit box, and an `ownership` list. By default the list is open to every agent at priority `0`. Tie one to an agent with:

```matlab
wp = waypoint('radius', 0.1, 'name', 'WP-alpha');
wp.SetGLOBAL('position', [10; 0; 0]);
wp.CreateAgentAssociation(agent, 5);   % higher number is preferred
```

`UpdateTargetWaypoint` keeps `this.targetWaypoint` on the visible waypoint with the highest priority that this agent has not already achieved. Achievement is local: the body-frame range is smaller than the sum of the two radii (`GetTargetCondition`). Achieved ids are stored on `this.achievedWaypoints`. When the visible set is empty, or every associated waypoint is achieved, the agent sets `GLOBAL.idleStatus` to `true`.

The simulator's own waypoint event is separate. It fires when the global collision test says the agent overlaps a waypoint it is allowed to claim. An agent can therefore mark a waypoint achieved in memory on the same step the event log records `eventType.waypoint`.

`GetTargetHeading` returns the unit vector from the body origin toward `targetWaypoint.position`. With no target, that heading is body `+x`.

## Obstacles

`obstacle` is a passive `objectDefinition`:

- `type` is `obstacle` and the hit box is spherical.
- It has no detection radius and receives no observation packet.
- `main` inherits the `objectDefinition` default, which integrates a zero input through `SingleIntegratorDynamics` and publishes it with `GlobalUpdate`. The published velocity is zero, so an initial scenario velocity is integrated for one step and the obstacle then holds position.

`obstacle_cuboid` and `obstacle_spheroid` replace the mesh. Leave `main` alone for a static obstacle. Override `main` and keep publishing a non-zero `GLOBAL.velocity` for an obstacle that moves on a scripted path.

## Sensors and dynamics containers

`agent.SENSORS` defaults to perfect measurements:

```matlab
range            = inf
sigma_position   = 0
sigma_velocity   = 0
sigma_rangeFinder = 0
sigma_camera     = 0
sampleFrequency  = inf
```

`sampleFrequency` is capped at the simulation frequency `1/dt` during setup. `agent.DYNAMICS` holds `maxLinearVelocity`, `maxLinearAcceleration`, `maxAngularVelocity`, and `maxAngularAcceleration`, sized `3x1` in 3D and reduced in 2D.

`agent.DATA` is a free struct for anything you want saved on the object and later read by the figure generator. `writeAgentData` is the convention used by the computation-time figure.

## Minimal agent

This subclass senses with the default packet, flies a constant body velocity, and publishes a 12-state update. It is the same contract as `agent.main`, with the decision block left open.

```matlab
classdef agent_constant < agent
    methods
        function [this] = agent_constant(varargin)
            this@agent(varargin);
            this.radius = 0.5;
            this.detectionRadius = 25;
            this = this.ApplyUserOverrides(varargin);
        end

        function [this] = main(this, ENV, varargin)
            [this, ~, ~] = this.GetAgentUpdate(ENV, varargin{1});

            % Body-frame rates: forward 0.5 m/s, yaw 0.1 rad/s.
            this.localState(7:9)  = [0.5; 0; 0];
            this.localState(10:12) = [0; 0; 0.1];

            this = this.GlobalUpdate_3DVelocities(ENV.dt, this.localState);
        end
    end
end
```

`objects/agents/agent_example.m` is the in-tree sketch of the same idea: a constructor, a `main` that reads `ENV` and the packet, a local integrator, and a global publish. The methods that exist on the base classes today are `GetAgentUpdate`, `GlobalUpdate`, `GlobalUpdate_3DVelocities`, and `GlobalUpdate_direct`. `examples/setup_example.m` shows the surrounding script: build a cell array of `agent_example`, place it with `GetScenario_concentricRing`, and pass the cell array to `OMAS_initialise`.

## Checklist for a new agent

- The file lives in `objects/agents/` and subclasses `agent` or `agent_2D`.
- The constructor calls the superclass, then `ApplyUserOverrides`.
- `setup` leaves `localState` consistent with `Is3D()` and records `priorState`.
- `main(this, ENV, varargin)` returns the updated object on every path.
- Observations are read from `varargin{1}` and passed to `GetAgentUpdate` when you want memory and waypoints.
- Some `GlobalUpdate*` method runs before `main` returns, so the velocity and quaternion it writes are integrated on the next step.
- `detectionRadius`, `radius`, and `v_max` are set for the scenario scale.
- `idleStatus` becomes `true` when the agent should let the simulation end.
