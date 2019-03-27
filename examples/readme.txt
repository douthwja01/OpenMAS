%%%% MATLAB MULTI-AGENT SIMULATOR (OpenMAS) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Author James A. Douthwaite (jadouthwaite1@sheffield.ac.uk)

%%%% GENERAL README %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OVERVIEW

This file is intended to provide an brief introduction to this open multi-agent
simulation (OpenMAS) tool. This directory contains a number of folders that 
together form a framework for the simulation of generic multi-agent scenarios.  
Inside this directory:
- data 
	+ The outputted simulation data and figures, ordered by simulation date 
	 and time of execution.

- environment 
	+ The directory containing the simulation functions and utilities.
		[it is advised you do not change anything within this folder]

- events    
	+ The simulations event set which are triggered throughout the runtime.
	  [it is advised you do not change anything within this folder]

- objects   
	+ The directory of user/simulation object definitions.

- scenarios 
	+ A directory of scenario definitions. Functions within here are 
		designed to generate the object initial conditions in the global
	  coordinate system. Scenario.fig, scenario.mat are auto generated.

IMPORTANT: Please ensure that all these directories are on the system path when 
	   creating a simulation setup function (see setup_example.m).

%%% SIMULATION PROCEDURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

An overivew of the complete process:
1. Add the set of directories to the path of your script [i.e 'addpath('objects')'] 
2. Define the matrix of initial states using the 'scenarioBuilder' or use and existing 
		scenario defined in /scenarios.
3. Initialise the array of objects/agents with these initial conditions.
4. Pass the array to the simulation with additional simulation parameters via the
   simulation wrapper function 'simulation_initialise'.
5. Find the output DATA and META structures, in addition to selected figures within 
   the output 'data' directory (labelled with the current date+time).

This process can be seen described in 'setup_example.m'.

%%% FIGURES AVAILABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The figures currently available for auto-generation are defined in 
'simuation_figureGenerator.m'. These figures will be generated from the data generated
in the current simulation. Providing the following arguements either individually or 
as a cell array of labels to 'simulation_initialise' will result in the following:

- 'ALL' 				 + Outputs all available figures.
- 'EVENTOVERVIEW + A summary of the simulation events that occured over the 
									 simulation time.
- 'COLLISIONS'   + Snapshots of the the collision instances between agents.
- 'TRAJECTORIES' + Produces a timeseries plot of each agents global state trajectories.
- 'SEPERATIONS'  + The inter-agent seperations from the perspective of each agent.
- 'INPUTS'			 + The control input trajectories, if stored to the agent.DATA property.
- 'ISOMETRIC'		 + An isometric summation of the global agent trajectories.
- '4VIEW'				 + A four perspective view of the scenario, annomated over time.
- 'GIF'					 + An animation of the isometric view.
- 'TIMES'				 + Plots of the computation times, if stored to the agent.DATA property.

example of passing figure requests to simulation_initialise:

[DATA,META] = simulation_initialise('objects',agentArray,'figures','ALL');

%%% MAIN FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The main function the user will interact with is the simulation wrapper function
named 'simulation_initialise.m' and your desired object class. The rest of the 
simulation is very much self-contained. Notes on some of the key functions are given
below:

- simulation_intialise 
	+ Handles the all the inputs to the simulation and returns 
		the output data.
	+ A full list of parameters can be found by inspecting the
		'validateSimulationInputs' function. 
	+ The main parameters include the object/agent and figure 
		request array passed as string pairs. For example:
		[DATA, META] = simulation_initialise('simtime',10,...
				       'objects',objectSet,'figures',listOfFigureLabels)]	 

- scenarioBuilder 
	+ Only needed if you wish to define a scenario not listed.
	+ A class object, initialised with a given number of objects. 
		[i.e testScenario = scenarioBuilder(numberofobjects)]
	+ Once initialised, the scenario builder methods can be called to
		design scenario description as a set of positions, velocities 
		and quaternion attitudes.
		[i.e scenarioConfig = testScenario.random()]
	+ Inspect this class to find the list of available scenarios.

- simulation_figureIndex 
	+ Contains the list of available output figures which can be
		selected via the 'simulation_initialise' function. Multiple
		figures can be requested as via an array of string names.
		[i.e. figureSet = {'eventoverview','trajectories'};]
	+ The simulation outputs the figures 'eventoverview' and 'fig' 
		by default.
	+ To request all figures, pass the input 'all' to figures.

%%% FINAL COMMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

To test and agent you define yourself, copy the example agent (agent_example) and 
begin building its functions into it. Build a cell array of agents and hand them to
a predefined scenario to get started and initialise their global properties. 

Finally define the simulation parameters (agents, time, timestep) and the requested
figures to be auto-generated.

Run this setup script and the output data will be returned the working directory.


		         