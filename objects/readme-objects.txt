%%%% OPENSOURCE MULTI-AGENT SIMULATOR (OpenMAS) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Author James A. Douthwaite (douthwaite@gmail.com)

%%% --> OBJECT CREATION AND DEFINTION 

OVERVIEW

The directory 'objects' contains a set of standard object definitions in the form 
of Matlab classes. Each object is considered iterations of the initial 
'objectDefinition' class, containing properties unique to that object or subsequent
object groups. Each object must provide certain methods for the simulation to interact
with, or rely on those of its superclasses, in order for the progress of the object to
be evaluated at the next timestep.

%%% KEY CLASS NAMES & HEIRARCHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%

- OBJECT DEFINITION (objectDefinition.m)
  + A matlab class containing numerous basic methods for interfacing with 
    the OpenMAS simulator (i.e global properties - global position, global velocity 
    and global quaternion pose).
  + Provides the incremental numbering/naming regime for unique identification.
  + Provides some generic tools updating global properties of objects with different 
  	state vectors. 
  + Provides basic object dynamic functions (particles, integrators, cars, bicycles)
  	to be imported.
  + Provides some basic vector manipulation utilities.
  + Provides generic geometry importation methods for .stl files by the same name as 
  	the class.

- AGENT (agent.m)
  + Inherits the 'objectDefinition' class.
  + Assumed to be 'non-passive', information on its surrounding is sent to the agent
  	each timestep.
  + The simulation then knows to look for the method named 'main', 
  	ENSURE THIS FUNCTION RETURNS THE UPDATED OBJECT.
  + Defines the a number of placeholder parameters such as 'detectionRange' which 
  	parameterise the agents visual range.

- 2D AGENT (agent_2D.m)
  + Inherits the 'agent' class.
  + Overrides the tools provided by 'agent.m' to facilitate simulation entirely in the
  	X/Y plane using 2D vectors.
  + The updateGlobalState function is important as it maps the resulting dynamics back
  	to the 3D environment.

OBSTACLE (obstacle.m)
  + Inherits the 'objectDefinition' class.
  + The base class for objects with a 'passive' main cycle. These objects
    will act passively, under constant conditions and will be assumed to move along a 
    given trajectory.

WAYPOINT (waypoint.m)
  + Inherits the 'objectDefinition' class.
  + A non-collidable object class. 
  + Waypoints can either be assigned an 'association' with a agent(s); meaning only 
  	they can observe the waypoint, otherwise it remains visible to all by default.
  + Used to define objectives in the environment; the way-point logic can be found in 
  	'agent.getAgentUpdate'.

OTHER
  + All other objects contained within the directory will be diviations upon the
    aforementioned object definitions.
  + Each object may in tern provide alternative procedures, algorithms for
    interpreting global or environmental information provided to it.
  + In addition to this, the user may provide a descrete dynamic model as part 
    of the 'main' method to observe non-linear or custom dynamical systems.

%%% CREATING NEW DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

There are several properties the user should be aware of in creating new 
simulation object definitions. An example of how create a basic agent can be 
found inside (agent_example.m) where this is shown more clearly, otherwise:

BRIEFLY:
- All objects
	+ Inherit the 'objectDefinition' class, properties and methods.
	+ Must compute their algorithm each cycle of the 'main' function.
	+ Define their 'obj.localState' using a defined or inherited dynamic function.
	+ Update their global properties in the .VIRTUAL structure using their new state
	  for the next timestep.

- All objects inheriting the 'agent' class must have a 'main' method
  which contains the evolution of the object over the timestep.
	+ Within this, the user should call a dynamic function such as the inherited 
	'dynamics_singleIntegrator' method or 'dynamics_simple' to define basic dynamics,
	 or define a custom dynamics function unique to that agent.
	+ Once a new state is defined, calling a 'updateGlobalProperties' 
	  method will update the objects global properties from its objects own
	  property 'localstate'.

- All objects must have the property VIRTUAL, which provides the static 
  properties for the simulator.
	+ .size - A characteristic dimension to register collisions within
		 (i.e. a radius defining a bubble around the object).
	+ colour - An RGB colour vector, unique or arbitrarily assigned.
	+ symbol - A unique symbol character. unique or arbitrarily assigned.
	+ geometry - A structure defining the geometry of the agent in local coordinates.

- All objects inheriting 'objectDefinition', can be initialised with a unique
  name string. Failing to provide one results in a string being generated, in
  addition to a unique objectID number.  

To begin creating an agent, start with the agent_example, observe the methods and 
properties it inherits and the options available before adding custom functions and 
properties which may already be defined. 
 


