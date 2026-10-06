# OpenMAS simulation guide

OpenMAS is a MATLAB framework for simulating multi-agent scenarios. A setup script builds a cell array of object handles, places them in a global scene, and hands that array to `OMAS_initialise`. The simulator then steps every object forward, records trajectories and events, and writes figures to a dated output folder.

The runnable walk-through is [`setup_example.m`](setup_example.m). How to write an agent is covered in [Agents and objects](agents.md).

## Repository layout

| Path | Role |
| --- | --- |
| `core/` | Simulator, analysis, figures, utils, and events. Treat this as the framework. |
| `core/utils/` | Shared MATLAB utilities used by the framework. |
| `objects/` | Nested object definitions (see below). |
| `objects/base/` | Base classes: `objectDefinition`, `agent`, `agent_2D`, `agent_3D`, `waypoint`. |
| `objects/agents/` | Concrete agent algorithms. |
| `objects/obstacles/` | Passive bodies (`obstacle`, planetoids, etc.). |
| `objects/vehicles/` | Vehicle models (`quadcopter`, `ARdrone`, …); legacy under `vehicles/legacy/`. |
| `objects/models/` | STL meshes named after their class. |
| `objects/tools/` | Shared helpers used by agent algorithms. |
| `scenarios/` | Functions that place an object set in the global frame, plus `scenario_builder`. |
| `examples/` | Setup scripts and these notes. |
| `studies/` | Longer study scripts. |
| `data/` | Default session output, one folder per run. |
| `toolboxes/` | Optional third-party libraries, including INTLAB. |

`core/OMAS_system.GetFileDependancies` adds the framework paths when a simulation starts. A setup script still needs `objects` (recursively) and `scenarios` on the path before it constructs classes:

```matlab
addpath('core');
addpath(genpath('objects'));
addpath('scenarios');
```

## Running a simulation

1. Construct each participant as a class handle (`agent_example`, `obstacle`, `waypoint`, and so on). Property overrides are name-value pairs, for example `agent_example('radius', 0.5, 'name', 'alpha')`.
2. Place those handles with a scenario function such as `scenario_concentric_ring`, or assign global pose yourself through `SetGLOBAL`.
3. Call `OMAS_initialise` with the object cell array and any timing or figure options.
4. Read the returned `DATA` and `META` structures. The same session is also written under the output path.

```matlab
agentIndex = cell(1, 3);
for k = 1:3
    agentIndex{k} = agent_example('radius', 0.5);
end

objectIndex = scenario_concentric_ring( ...
    'agents', agentIndex, ...
    'agentOrbit', 5, ...
    'agentVelocity', 2, ...
    'waypointOrbit', 10, ...
    'plot', true);

[DATA, META] = OMAS_initialise( ...
    'objects', objectIndex, ...
    'duration', 10, ...
    'dt', 0.1, ...
    'figures', {'events', 'plan', 'gif'}, ...
    'warningDistance', 2, ...
    'verbosity', 1);
```

`OMAS_initialise` closes open figures, checks the object vector, builds the `META` record, runs `OMAS_process`, then runs `OMAS_analysis`.

## The `OMAS_initialise` interface

Inputs are name-value pairs. Unrecognised names are ignored. Timing fields are applied to the internal `TIME` structure.

| Name | Default | Meaning |
| --- | --- | --- |
| `objects` | none | Cell array of object handles. Required, and every entry must be an object. |
| `duration` | `1000` | Simulation length in seconds, measured from `startTime`. |
| `dt` | `0.5` | Fixed step in seconds. Must be smaller than `duration`. |
| `startTime` | `0` | First sample time. The time vector is `startTime:dt:endTime`. |
| `idleTimeOut` | `1` | Extra seconds to keep running after every agent reports idle. |
| `figures` | `{'events','fig'}` | Figure labels, a string, or a cell array of strings. See below. |
| `warningDistance` | `5` | Proximity margin in metres, measured outside the two object radii. |
| `conditionTolerance` | `1e-3` | Numeric tolerance used by collision and detection tests. |
| `outputPath` | `<pwd>/data` | Directory that receives a new session folder. |
| `verbosity` | `2` | `0` is quiet, `1` prints each step, `2` also prints each object update. |
| `threadPool` | `false` | When `true`, agent cycles run in a `parfor` over a local pool. Forced off in Monte-Carlo mode. |
| `monteCarloMode` | `false` | Writes each cycle into a uniquely named folder and disables the thread pool. |
| `gui` | `true` | Logical flag reserved for interface use. |
| `publishMode` | `false` | Logical flag reserved for figure publishing. |

`duration` must be at least `dt`, `startTime` must be at or before the end time, and `idleTimeOut` must be zero or positive.

Global pose on every object is checked before the loop starts. `position` and `velocity` are `3x1` columns, and `quaternion` is a `4x1` column. The simulator then calls `setup(localVelocity, localEuler)` on each object, with the velocity and ZYX Euler angles expressed in that object's body frame.

Outputs:

- `DATA` holds the sampled trajectories, timing, and the event statistics added by analysis.
- `META` is the terminal simulation record: time, object summaries, and the output path.

## What one time step does

`OMAS_process` walks `TIME.timeVector`. At step `k` the order is:

1. **Idle exit.** If every agent has `idleStatus` set, the loop prints a countdown and stops once `idleTimeOut` has elapsed. A scene with no agents runs to `duration`.
2. **Stamp the clock.** `currentTime` and `currentStep` are updated.
3. **Record the state.** Each object's global state is copied into `DATA.globalTrajectories` at this step, before the objects move.
4. **Integrate the world model.** `OMAS_updateGlobalStates` advances every meta object with the velocity and quaternion currently stored on that object:

   `position = position + velocity * dt`

   Attitude `R` is rebuilt from the quaternion. The position written inside the object is not copied across. Motion in the scene is whatever velocity `main` published on the previous step.
5. **Separations and events.** Pairwise centroid separations, warning margins, collisions, and detections are evaluated from that integrated state. New and cleared conditions become `eventDefinition` records.
6. **Object cycles.** Each object is updated. Agents receive an observation packet. Obstacles and waypoints are stepped with time only. Objects of type `misc` are left as they are. Velocity and quaternion written here are integrated on the next step.
7. **Advance.** The step index increments. The first recorded column is the initial condition. Later columns are the state integrated from the velocity published so far.

Global coordinates are East-North-Up, stored as MATLAB `x`, `y`, `z`. The rotation matrix `R` on each meta object maps a global vector into that object's body frame (`body = R * global`). Body `x` is forward.

### Observations

Only objects whose `GLOBAL.type` is `OMAS_objectType.agent` receive observations. For every other object that is currently in detection, the simulator builds one struct and concatenates them into the packet passed to `main`:

| Field | Contents |
| --- | --- |
| `objectID`, `name`, `type`, `colour` | Identity copied from the observed object. |
| `position`, `velocity` | Relative pose in the **observer's body frame**. 3-vectors for a 3D agent, 2-vectors for an `agent_2D`. |
| `range`, `heading`, `elevation`, `width` | Spherical range, azimuth, elevation, and apparent angular width. |
| `radius` | Characteristic radius used for the angular-width estimate. |
| `geometry` | Vertices, faces, normals, and centroid, rotated into the observer frame. Empty when the object has no mesh. |
| `DEBUG` | `globalPosition`, `globalVelocity`, and waypoint `priority`. This is ground truth for debugging. |

`position` is `p_other - p_agent`, rotated by the observer's `R`. An empty packet means nothing is inside the detection volume this step.

Detection itself is decided before `main` runs:

- The observer's `detectionRadius` sphere is tested against the other object's sphere (`radius`).
- A waypoint is reported only when `waypoint.IDAssociationCheck` accepts the observer.
- If the observed object has a triangular mesh, a face must intersect the detection sphere before the detection sticks.
- Detection is one-way. Agent A can see agent B while B's shorter horizon misses A.

The first time a pair enters or leaves detection, warning, collision, or waypoint range, an event is logged. Repeating the same condition on later steps does not emit another event until the condition clears (`null_detection`, `null_warning`, `null_collision`).

Warnings use the gap between centroids after both radii are subtracted, compared with `warningDistance`. Collisions use `OMAS_collisionDetection` and the object's hit box (`OMAS_hitBoxType`: `none`, `spherical`, `AABB`, `OBB`, `capsule`, or `mesh`). Waypoint capture is a collision test against a waypoint the agent is allowed to claim, and it does not also raise a collision event.

### The `ENV` structure passed into `main`

`UpdateObjects` passes `ENV = META.TIME` plus `ENV.outputPath`. The fields an agent can read are:

| Field | Meaning |
| --- | --- |
| `dt` | Step length, seconds. |
| `currentTime` | Time of the sample just recorded. |
| `currentStep` | Index into `timeVector`, starting at 1. |
| `timeVector` | Full planned time row. |
| `startTime`, `endTime`, `duration` | Schedule bounds. `endTime` is rewritten to the actual stop time when the loop finishes. |
| `numSteps`, `frequency` | Step count and `1/dt`. |
| `idleTimeOut` | Idle hold, seconds. |
| `outputPath` | Session directory. |

Agents are called as:

```matlab
this = this.main(ENV, observationPacket);
```

Obstacles and waypoints are called as `this.main(ENV)` with no packet. The method must return the same object. On the following step the simulator reads `GetGLOBAL('velocity')`, `GetGLOBAL('quaternion')`, and `GetGLOBAL('idleStatus')`, and moves the body by `velocity * dt`. Editing `localState` alone leaves the world state unchanged. See [Agents and objects](agents.md).

## Figures

Pass labels to the `figures` argument. The match is case-insensitive. `'all'` expands to the full set. `'none'` suppresses output. Any other string is skipped with a warning.

| Label | Output |
| --- | --- |
| `events` | Event timeline (detections, warnings, collisions, waypoints). |
| `collisions` | Snapshots of collision instants. |
| `trajectories` | Time series of each object's global state. |
| `separations` | Inter-object separation from each object's point of view. |
| `closest` | Minimum separation over the run. |
| `avoidance` | Avoidance-rate summaries. |
| `inputs` | Control inputs stored on `agent.DATA`. |
| `plan` | Top-down plan view of the trajectories. |
| `isometric` | 3D isometric trajectories. |
| `gif` | Animated isometric GIF. |
| `avi` | Animated isometric AVI. |
| `times` | Per-step computation times stored on `agent.DATA`. |

`agent.writeAgentData(TIME, loopIndicator, loopDuration)` is the helper that fills `DATA.indicator`, `DATA.dt`, `DATA.steps`, and `DATA.time` for the `inputs` and `times` figures.

## Files written for each session

`ConfigureOutputDirectory` creates `<outputPath>/<yyyy-mm-dd@HH-MM-SS>-session-data/`. Analysis then saves:

| File | Contents |
| --- | --- |
| `META.mat` | Simulation configuration and per-object status at the end of the run. |
| `DATA.mat` | `timeVector`, `dt`, `stateIndex`, `globalTrajectories`, and event statistics. |
| `EVENTS.mat` | Vector of `eventDefinition` objects. |
| `OBJECTS.mat` | The final object handles, including agent `MEMORY` and `DATA`. |

`DATA.globalTrajectories` has one column per planned step and `10 * nObjects` rows. For object `i`, rows `stateIndex(1,i):stateIndex(2,i)` are `[x; y; z; vx; vy; vz; qw; qx; qy; qz]`. Early termination leaves trailing `NaN` columns. The position in that column is the simulator integral of published velocity, sampled at the start of the step.

`META.OBJECTS(i)` is the simulator's public record of object `i`: `objectID`, `name`, `class`, `type`, `hitBox`, `colour`, `symbol`, `radius`, `detectionRadius`, `idleStatus`, `globalState` (`[position; velocity; quaternion]`), `R`, `relativePositions`, and `objectStatus`.

## Scenarios

A scenario function accepts the agent cell array, writes each object's global position, velocity, and quaternion, and usually appends waypoints. `scenario_concentric_ring` is the pattern used by `setup_example.m`: agents on a ring, opposing waypoints, optional planar position noise.

`scenario_builder` generates the underlying point sets in the same East-North-Up frame. Construct it with the object count, then call a generator:

| Method | Layout |
| --- | --- |
| `planarDisk` | Concentric rings in a plane. |
| `planarRing` | A single ring. |
| `planarAngle` | A planar angular sector. |
| `regularSphere` | Nodes on a sphere. |
| `helix` | A helical curve. |
| `line` | Evenly spaced along a line. |
| `random`, `randomNormal`, `randomUniform`, `randomSphere` | Random placements. |

Each generator fills `.positions` (`3 x n`), `.velocities` (`3 x n`), and `.quaternions` (`4 x n`). Copy those columns onto the objects with `SetGLOBAL`, then return one cell array of every participant. `scenario_builder.plotObjectIndex(objectIndex)` draws the initial scene.

Ready-made wrappers live in `scenarios/scenario_*.m` (concentric ring and sphere, two lines, corridor, waypoint curve, obstacle track, formation splits, random fields, and the Earth-orbit study).

## Monte-Carlo studies

`OMAS_monteCarlo` repeats many OpenMAS cycles. A **study** is a set of sessions, a **session** is one object population (typically one algorithm), and a **cycle** is one perturbed run of that population.

```matlab
mc = OMAS_monteCarlo( ...
    'objects', studyObjects, ...   % cell array, sessions by objects
    'cycles', 5, ...
    'isParallel', false, ...
    'directory', outputPath);
mc.EvaluateAllCycles();
```

The object matrix is `sessions x objects`. Each cycle perturbs global positions before calling the same simulator. `examples/setup_monte_carlo.m` builds that matrix for several avoidance algorithms and population sizes.

## Object types

`OMAS_objectType` is the label stored on `GLOBAL.type`:

| Enum | Value | Stepped by |
| --- | --- | --- |
| `misc` | 0 | Nothing. The object stays at its last global state. |
| `agent` | 1 | `main(ENV, observationPacket)` |
| `obstacle` | 2 | `main(ENV)` |
| `waypoint` | 3 | `main(ENV)` |

`eventType` labels the rows of the event log: `detection`, `warning`, `collision`, `waypoint`, and the matching `null_*` clear events.

## Object dictionary

Classes that can be placed in a scenario. `Inherits` is the parent list from the `classdef` line. Paths are relative to `objects/`. The methods an agent must implement are described in [Agents and objects](agents.md).

### Bases

| Class | Inherits | File | Description |
| --- | --- | --- | --- |
| `agent_tools` | `objectDefinition` | `base/agent_tools.m` | Shared memory, waypoint selection, and global-pose updates. Mixed into `agent`. |
| `agent` | `objectDefinition`, `agent_tools` | `base/agent.m` | Active agent base. Sensing, a detection radius, kinematic limits, and `main(ENV, observations)`. |
| `agent_2D` | `agent` | `base/agent_2D.m` | Planar agent. Six-state local vector `[x; y; psi; dx; dy; dpsi]` and 2D observations. |
| `agent_3D` | `agent` | `base/agent_3D.m` | Spatial kinematic agent. Twelve-state local vector and the simple and PID velocity trackers. |

### Templates

| Class | Inherits | File | Description |
| --- | --- | --- | --- |
| `agent_example` | `agent` | `agents/agent_example.m` | Skeleton agent. Constant body rates, a local integrator, and a global pose update. |
| `agent_2D_test` | `agent_2D` | `agents/agent_2D_test.m` | Planar sandbox with a shortened local state for trying 2D updates. |

### Velocity obstacles (3D)

| Class | Inherits | File | Description |
| --- | --- | --- | --- |
| `agent_VO` | `agent` | `agents/agent_VO.m` | Geometric velocity-obstacle avoidance. Samples a feasible velocity grid and selects a course outside neighbouring velocity obstacles. |
| `agent_RVO` | `agent_VO` | `agents/agent_RVO.m` | Reciprocal velocity obstacles (van den Berg et al., 2008). Each agent takes half of the pairwise avoidance. |
| `agent_HRVO` | `agent_RVO` | `agents/agent_HRVO.m` | Hybrid reciprocal velocity obstacles. Chooses the VO or RVO cone from which side of the obstacle the current velocity lies on. |

### Velocity obstacles (2D)

| Class | Inherits | File | Description |
| --- | --- | --- | --- |
| `agent_2D_VO` | `agent_2D`, `agent_VO` | `agents/agent_2D_VO.m` | Planar velocity-obstacle avoidance. |
| `agent_2D_RVO` | `agent_2D_VO`, `agent_RVO` | `agents/agent_2D_RVO.m` | Planar reciprocal velocity obstacles. |
| `agent_2D_HRVO` | `agent_2D_RVO`, `agent_HRVO` | `agents/agent_2D_HRVO.m` | Planar hybrid reciprocal velocity obstacles. |
| `agent_2D_ORCA` | `agent_2D_RVO` | `agents/agent_2D_ORCA.m` | Optimal reciprocal collision avoidance. Builds ORCA half-plane constraints over a finite time horizon. |

### Vector sharing and intervals

| Class | Inherits | File | Description |
| --- | --- | --- | --- |
| `agent_vectorSharing` | `agent` | `agents/agent_vectorSharing.m` | Geometric avoidance that shares the separating correction between the two agents in a conflict. |
| `agent_2D_vectorSharing` | `agent_2D`, `agent_vectorSharing` | `agents/agent_2D_vectorSharing.m` | Planar vector-sharing avoidance. |
| `agent_interval` | `agent` | `agents/agent_interval.m` | Interval-analysis helpers (INTLAB) mixed into the interval avoidance agents. |
| `agent_IA` | `agent_vectorSharing`, `agent_interval` | `agents/agent_IA.m` | 3D interval avoidance. Vector sharing with interval bounds on the uncertain relative state. |
| `agent_2D_IA` | `agent_2D_vectorSharing`, `agent_interval` | `agents/agent_2D_IA.m` | Planar interval avoidance. |

### Formation

| Class | Inherits | File | Description |
| --- | --- | --- | --- |
| `agent_formation` | `agent` | `agents/agent_formation.m` | Formation tracking from an adjacency matrix. Supplies the bearing and displacement control laws used by the formation variants. |
| `agent_formation_boids` | `agent_formation` | `agents/agent_formation_boids.m` | Reynolds boids. Separation, alignment, cohesion, and waypoint seeking, weighted and sent to the kinematic controller. |
| `agent_formation_VO` | `agent_VO`, `agent_formation` | `agents/agent_formation_VO.m` | Velocity-obstacle avoidance while tracking a formation. |
| `agent_formation_RVO` | `agent_formation_VO`, `agent_RVO` | `agents/agent_formation_RVO.m` | Reciprocal velocity-obstacle avoidance while tracking a formation. |
| `agent_formation_HRVO` | `agent_formation_RVO`, `agent_HRVO` | `agents/agent_formation_HRVO.m` | Hybrid reciprocal avoidance while tracking a formation. |
| `agent_2D_formation_VO` | `agent_2D_VO`, `agent_formation` | `agents/agent_2D_formation_VO.m` | Planar velocity-obstacle formation agent. |
| `agent_2D_formation_RVO` | `agent_2D_formation_VO`, `agent_2D_RVO` | `agents/agent_2D_formation_RVO.m` | Planar reciprocal-velocity-obstacle formation agent. |
| `agent_2D_formation_HRVO` | `agent_2D_formation_RVO`, `agent_2D_HRVO` | `agents/agent_2D_formation_HRVO.m` | Planar hybrid-reciprocal formation agent. |
| `agent_2D_formation_ORCA` | `agent_2D_formation_RVO`, `agent_2D_ORCA` | `agents/agent_2D_formation_ORCA.m` | Planar ORCA avoidance while tracking a formation. |

### Vehicles

| Class | Inherits | File | Description |
| --- | --- | --- | --- |
| `quadcopter` | `agent` | `vehicles/quadcopter.m` | Quadrotor. Rotor plant in `DYNAMICS` and a position and velocity controller. |
| `quadcopter_legacy` | `quadcopter` | `vehicles/legacy/quadcopter_legacy.m` | Earlier quadrotor from Shiyu Zhao, mapped onto the quadcopter class with its original plant left intact. |
| `quadcopter_formation` | `quadcopter`, `agent_formation` | `vehicles/quadcopter_formation.m` | Quadrotor plant that also tracks an adjacency-matrix formation. |
| `ARdrone` | `quadcopter` | `vehicles/ARdrone.m` | Parrot AR.Drone. Outdoor configuration, battery ratings, sensors, and a body mesh on the quadrotor plant. |
| `ARdrone_new` | `quadcopter` | `vehicles/ARdrone_new.m` | Later AR.Drone parameter set. Indoor or outdoor configuration and battery ratings on the quadrotor plant. |
| `ARdrone_prev` | `quadcopter` | `vehicles/legacy/ARdrone_prev.m` | Earlier AR.Drone dynamic model kept for the LQR and MPC controllers. |
| `ARdrone_LQR` | `ARdrone_prev` | `vehicles/ARdrone_LQR.m` | Linear-quadratic regulator for velocity and position on the earlier AR.Drone plant. |
| `ARdrone_MPC` | `ARdrone_LQR` | `vehicles/ARdrone_MPC.m` | Model-predictive controller over a short horizon, built on the LQR drone. |
| `ARdrone_formation` | `ARdrone_LQR`, `agent_formation` | `vehicles/ARdrone_formation.m` | LQR AR.Drone that also tracks an adjacency-matrix formation. |
| `fixedWing` | `agent` | `vehicles/fixedWing.m` | Placeholder fixed-wing dynamics. Dubins-style state `[x; y; z; v; psi; theta]`. |
| `A10` | `fixedWing` | `vehicles/A10.m` | A-10 airframe on the fixed-wing plant. Geometry follows the class name. |
| `boeing737` | `fixedWing` | `vehicles/boeing737.m` | Boeing 737 airframe on the fixed-wing plant. Geometry follows the class name. |
| `globalHawk` | `fixedWing` | `vehicles/globalHawk.m` | Global Hawk airframe on the fixed-wing plant. Geometry follows the class name. |
| `ISS` | `agent` | `vehicles/ISS.m` | International Space Station as an orbiting body, with an orbital radius and orbital speed. |

### Obstacles

Passive bodies. The simulator steps them with `main(ENV)` and no observation packet. The default cycle publishes zero velocity, so an initial velocity is integrated for one step and the obstacle then holds position.

| Class | Inherits | File | Description |
| --- | --- | --- | --- |
| `obstacle` | `objectDefinition` | `obstacles/obstacle.m` | Generic obstacle. Spherical hit box and a six-state local vector. Static unless `main` is overridden to keep publishing a velocity. |
| `obstacle_spheroid` | `obstacle` | `obstacles/obstacle_spheroid.m` | Sphere mesh built from the radius when no geometry file is loaded. Spherical hit box. |
| `obstacle_cuboid` | `obstacle` | `obstacles/obstacle_cuboid.m` | Axis-scaled cuboid mesh with an oriented bounding-box hit box. `Xscale`, `Yscale`, and `Zscale` set the side lengths. |
| `planetoid` | `obstacle_spheroid` | `obstacles/planetoid.m` | Spinning spheroid for satellite scenes. Stores inclination, orbit, orbital speed, mass, and a constant axial spin rate. |
| `earth` | `planetoid` | `obstacles/earth.m` | Earth parameters on the planetoid: mean radius 6.371e6 m and sidereal spin 7.2921150e-5 rad/s. |
| `moon` | `planetoid` | `obstacles/moon.m` | Earth's moon. Radius, mass, orbital radius, orbital speed, and a 5.145 degree inclination. |

### Waypoints

| Class | Inherits | File | Description |
| --- | --- | --- | --- |
| `waypoint` | `objectDefinition` | `base/waypoint.m` | Spherical goal. An ownership list ties it to one or more agents with a priority. Unowned waypoints are open to every agent at priority 0. Reaching one logs a waypoint event for an associated agent. |
