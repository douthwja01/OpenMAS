%%%% MATLAB MULTI-AGENT SIMULATOR (MAS) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Author James A. Douthwaite (jadouthwaite1@sheffield.ac.uk)

%%% --> OBJECT CREATION AND DEFINTION 

OVERVIEW

The directory 'objects' contains the definitions of various simulation classes. 
This is done in the form of matlab classes, containing each algorithm and physical
properties. Each object must provide certain methods for the simulation to interact
with in order for it to evaluate the progress of the agent at the next timestep. .

%%% KEY CLASS NAMES & HEIRARCHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%

- OBJECT DEFINITION (objectDefinition.m)
  + A matlab class containing numerous basic methods for integration with 
    the simulator.
  + Provides the incremental numbering/naming regime for unique identification.
  + Provides methods for initialising objects with given global properties (i.e 
    global position, global velocity and quaternion).
  + Provide a method for defining an objects local state.
  + Provides a basic particle motion model, inherited by all other objects by 
    default.
  + Provides some basic vector manipulation tools.

- AGENT (agent.m)
  + Inherits the 'objectDefinition' class.
  + Provides a basis for all objects with active 'cycles' (i.e are non-passive). 
  + The simulation then knows to look for the method named
    'processTimeCycle', ENSURE THIS FUNCTION RETURNS THE UPDATED OBJECT.
  + Defines the placeholder variables 'sensorRange' and 'sampleFrequency', 
    emulating an agents limited visual range and rate of computation.

OBSTACLE (obstacle.m)
  + Inherits the 'objectDefinition' class.
  + The base class for objects without a 'process time-cycle'. These objects
    will act passively and will be assumed to move along a given trajectory.

OTHER
  + All other objects contained within the directory will be diviations upon the
    aforementioned classes.
  + Each object may in tern provide alternative procedures, algorithms for
    interpreting global information provided to it via the simulation.
  + In addition to this, the user may provide a descrete dynamic model as part 
    of the 'processTimeCycle' method to observe interations on dynamic systems.

%%% CREATING NEW DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

There are several properties the user should be aware of in creating new 
simulation object definitions. An example of how create a basic agent can be 
found inside (agent_example.m) where this is shown more clearly, otherwise:

- All objects inheriting the 'agent' class must have a 'processTimeCycle' method
  which contains the evolution of the object over the timestep.
	+ Within this, the user should call the inherited 'stateDynamics_acceleration' method
	  or 'stateDynamics_velocities' to apply linear dynamics, or overwrite this 
    method with a process unique to that agent.
	+ Once a new state is defined, calling the 'updateGlobalProperties' 
	  method will update the objects global properties from its objects own
	  property 'localstate'.

- All objects must have the property VIRTUAL, which provides the static 
  properties for the simulator.
	+ .size - A characteristic dimension to register collisions within
		 (i.e. a radius defining a bubble around the object).
	+ colour - An RGB colour vector, unique or arbitrarily assigned.
	+ symbol - A unique symbol character. unique or arbitrarily assigned.

- All objects inheriting 'objectDefinition', can be initialised with a unique
  name string. Failing to provide one results in a string being generated, in
  addition to a unique objectID number.  
 


