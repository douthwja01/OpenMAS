%% THE OPENMAS TIMECYCLE (OMAS_mainCycle.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the progressive cycles of the simulation, for the
% environment, obstacles, agents and objects withing the simulation space.

% Author: James A. Douthwaite 06/10/2016

function [DATA,SIM,EVENTS,objectIndex] = OMAS_mainCycle(META,objectIndex)
% INPUTS:
% objectIndex - The cell array of object classes
% OUTPUTS:
% DATA        - The output DATA structure with timeseries data
% SIM         - The terminal META structure copy

global plotnum

%% PRE-SIMULATION DECLARATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
META.phase = 'META'; % Declare the simuation active phase
EVENTS = [];         % Reset the EVENT log container
% DETERMINE PLOT PROPERTIES
if ~exist('plotnum','var') || isempty(plotnum)
    plotnum = 1;    % Default to first plot
end
% PREPARE THE OUTPUT DATA CONTAINER
[DATA] = getOutputStructure(META); 

%% BEGIN TIME-STEP INTERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n[%s]\tLAUNCHING SIMULATION...\n\n',META.phase);
tic;
[META,objectIndex,DATA,EVENTS] = processTimeSeries(META,objectIndex,DATA,EVENTS);
DATA.computationTime = toc; % Get total elapsed simulation time
fprintf('\n[%s]\tSIMULATION COMPLETE (time elapsed: %ss)\n',...
        META.phase,num2str(DATA.computationTime));

% IF SOMETHING HAS BEEN PLOTTED, ITERATE PLOT NUMBER
if ~isempty(findall(0,'Type','Figure'))
    display('PLOTS OPEN');
    plotnum = plotnum + 1; % Increment to avoid overwriting
end    

%% FINALISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
META.phase = 'OUTPUT';                                                     % Change to OUTPUT phase
clearvars -except META EVENTS DATA objectIndex                             % Clear loose simulation variables
fprintf('[%s]\tDumping initial simulation variables to output directory...\n',META.phase);
sendOutputToFiles(META,EVENTS,DATA,objectIndex);
fprintf('[%s]\tClosing simulation...\n',META.phase);
SIM = META;
end

%% MAIN CYCLE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS TIME VECTOR
function [META,objectIndex,DATA,EVENTS] = processTimeSeries(META,objectIndex,DATA,EVENTS)
% This function computes the simulation loop cycle across the given
% descrete META.timevector:
% INPUTS:
% META - The initialised META data-structure
% objectIndex - The cell vector of objects (and/or agents)
% DATA        - The initialised output DATA structure
% EVENTS      - The initialised EVENTS structure list

% PROCESS TIME CYCLE
for step = 1:META.TIME.numSteps
    %% 1. ////////// UPDATE TIMESTEP PARAMETERS ///////////////////////////
    META.TIME.currentTime = META.TIME.timeVector(step);                    % The current sim-time
    META.TIME.currentStep = step;                                          % The current sim-step
    if META.verboseMode
        fprintf('[%s]\tStep: %s\tTime: %ss\n',META.phase,num2str(META.TIME.currentStep),num2str(META.TIME.currentTime));
    end
    
    %% 2. ////////// COLLECT GLOBAL AGENT/OBSTACLE STATES /////////////////
    for ID1 = 1:META.totalObjects
        % Collect the META.OBJECTS.state data (as the simulations understanding of state) 
        % and save in the output DATA.globalTrajectories set, must be done in
        % synchronously).
        DATA.globalTrajectories(DATA.stateIndex(1,ID1):DATA.stateIndex(2,ID1),META.TIME.currentStep) = META.OBJECTS(ID1).globalState;
    end
    
    %% 3, ///////// UPDATE VISUAL REPRESENTATIONS /////////////////////////
    if META.visualisation
        for ID1 = 1:META.totalObjects
            % PASS THE GLOBAL OBJECT DATA TO THE VR SCENE UPDATER
            OMAS_visuals(META.OBJECTS(ID1));
        end
    end
    
    %% 4. ////////// UPDATE SIMULATION/ENVIRONMENTAL META DATA ////////////
    objectSnapshot = objectIndex;                                          % Make a temporary record of the object set
    for ID1 = 1:META.totalObjects                                          % For each META object element
        % UPDATE EACH META.OBJECT AGAINST THE CURRENT META AND OBJECT
        % INDEX SNAPSHOT
        [META.OBJECTS(ID1),metaEVENTS] = updateSimulationMeta(META,objectSnapshot,objectIndex{ID1});  % Update META snapshot with equivilant objectIndex.state
        % LOG THE META EVENTS
        if ~isempty(metaEVENTS)                                            % If META events occurred this timestep
            EVENTS = vertcat(EVENTS,metaEVENTS);                           % Append to global EVENTS
        end
    end
    
    %% 5. ////////// COMPUTE AGENT CYCLES /////////////////////////////////
    if META.threadPool ~= 0
        objectSnapshot = objectIndex;                                      % Make a temporary record of the object set
        parfor (ID1 = 1:META.totalObjects)
            % MOVE THROUGH OBJECT INDEX AND UPDATE EACH AGENT
            [objectIndex{ID1},objectEVENTS] = updateObjects(META,objectSnapshot,objectIndex{ID1});            
            % LOG THE OBJECT EVENTS
            if ~isempty(objectEVENTS)
                EVENTS = vertcat(EVENTS,objectEVENTS);
            end
        end
    else
        for ID1 = 1:META.totalObjects
            % MOVE THROUGH OBJECT INDEX AND UPDATE EACH AGENT
            [objectIndex{ID1},objectEVENTS] = updateObjects(META,objectIndex,objectIndex{ID1}); % Update objectIndex snapshot with new META data
            % LOG THE OBJECT EVENTS
            if ~isempty(objectEVENTS)                                      % If objectEVENTS occur in this timestep
                EVENTS = vertcat(EVENTS,objectEVENTS);                     % Append to global EVENTS
            end
        end
    end
end
end

%% SIMULATION SUB-OPERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE THE SIMULATION META VARIABLES FROM OBJECT CHANGES
function [METAObjUpdate,metaEVENTS]     = updateSimulationMeta(SIM,objectProofSet,entity)
% This function updates the SIM.OBJECT given the new entity(object) state.
% The function moves through the current SIM.OBJECT set and update the
% simulations knowledge of that agent. Proximities are also evaluated for collision
% before the next agent cycle computations. Returns the updated META object
% and the list of events that have occurred.
% INPUTS:
% SIM            - A snapshot of the META structure
% objectProofSet - A snapshot of the objectIndex (access ownership properties)
% entity         - The object index element
%   .objectID
%   .quaternion
%   .globalVelocity
% OUTPUTS:
% newMETAobj     - The updated META.OBJECTS(i) object
% metaEVENTS     - The EVENT vector of META events

% DECLARATIONS
metaEVENTS = [];                                                           % Reset meta Events container
METAObjUpdate  = SIM.OBJECTS(SIM.globalIDvector == entity.objectID);       % Get the current META object associated with 'entity'
detectionRange = SIM.visabilityFactor*(METAObjUpdate.detectionRange);      % Get the detection range from newMETAObj

%% REALLOCATE THE GLOBAL PROPERTIES FOR 'ENTITY' AND SAVE TO META.OBJECT
% UPDATE THE GLOBAL TO BODY & BODY TO GLOBAL ROTATION MATRIX FROM THE NEW 
% GLOBAL TO ROTATED BODY QUATERNION
[METAObjUpdate.R_BG,METAObjUpdate.R_GB] = OMAS_axisTools.quaternionToRotationMatrix(entity.quaternion);

% REDEFINE GLOBAL DESCRIPTION
% The rotation matrix is updated each timestep and from the
% globalToBodyQuaternion.
% ENUlocalVelocity = OMAS_axisTools.NEDtoENU(entity.localState(4:6));
globalVelocity = entity.globalVelocity;                                       % TAKE THE GLOBAL VELOCITY DIRECTLY
globalPosition = METAObjUpdate.globalState(1:3) + globalVelocity*SIM.TIME.dt; % Calculate the new global position
globalToRotatedBodyQuaternion = entity.quaternion;                            % Get the global triad to rotated body rotation

% REBUILD GLOBAL STATES (ENTITY & META)
METAObjUpdate.globalState = [globalPosition;globalVelocity;globalToRotatedBodyQuaternion];

%% ITERATE THROUGH THE SIMULATION OBJECT SET
for index = 1:SIM.totalObjects                                                % Loop each agent record in the META structure
    % UPDATE GLOBAL PROPERTIES (META DATA)
    % GET THE OBJECT TO EVALUATE AGAINST "METAObjUpdate" 
    METAObjEvaluation = SIM.OBJECTS(index);
    % NEGLECT SELF EVALUATION
    if METAObjEvaluation.objectID == METAObjUpdate.objectID                   % Where not own ID1
       continue 
    end

    % THE NEW GLOBAL STATE (METAObjUpdate) IS NOW RE-EVALUATED AGAINST THE OTHER SIM.OBJECT DATA
    % Define the new relative global state (x_rel = x_b - x_a)
    relativePosition = METAObjEvaluation.globalState(1:3) - METAObjUpdate.globalState(1:3);          % Relative global position
    METAObjUpdate.relativePositions(index,:) = relativePosition;                                     % Allocate new obj relative position vector
    METAObjUpdate.euclideanSeperation(index,1) = abs(sqrt(sum(relativePosition.^2)));                % Scalar seperation value
    euclideanSeperation = METAObjUpdate.euclideanSeperation(index,1);                                % Get the seperation for validation
    geometricLimits = METAObjUpdate.size + METAObjEvaluation.size;                                   % Evaluate the combined agent radii
    
    % CHECK IF ENTITY IS WAYPOINT; CANNOT CREATE META EVENTS, MOVE TO NEXT OBJECT
    if METAObjUpdate.type == OMAS_objectType.waypoint
        SIM.OBJECTS(SIM.globalIDvector == METAObjUpdate.objectID) = METAObjUpdate;                   % Assign the updated SIM.OBJECT
        continue 
    end
    
    % CHECK IF THE AGENT IS BEING UPDATED AGAINST IS A WAYPOINT'S META
    if METAObjEvaluation.type == OMAS_objectType.waypoint
        % DETERMINE IF THE REFERENCE META objectID APPEARS IN THE META EVAL
        waypointOwnership = [objectProofSet{index}.ownership.objectID];                              % Get the ID set for the entity
        % AUTHORISE DETECTION IS VALID OWNER OR HAS NO OWNERS
        authorisedDetection = any(ismember(waypointOwnership,METAObjUpdate.objectID)) || any(isnan(waypointOwnership));
        isCollidable = 0;
    else
        % THE AGENT IS BEING UPDATED AGAINST ANOTHER AGENT
        authorisedDetection = 1;
        isCollidable = 1;
    end
    
    %% ASSESS DETECTIONS CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detections must be made on the basis of visual range/proximity and 
    % authorisation in order maintain generality.
    objectStatus = METAObjUpdate.objectStatus(index,eventType.detection);                          % Has the object already been detected
    detectionCondition = detectionRange >= euclideanSeperation && objectStatus == 0;               % Detection condition (all objects)
    detectionNullCondition = detectionRange < euclideanSeperation && objectStatus == 1;            % Detection loss condition (all objects)
    if detectionCondition && authorisedDetection
        % GENERATE DETECTION EVENT
        % Create EVENT between the updated and problem META objects
        [METAObjUpdate,EVENT] = OMAS_eventHandler(SIM,METAObjUpdate,METAObjEvaluation,eventType.detection);
        metaEVENTS = vertcat(metaEVENTS,EVENT);                                                    % Add new META EVENTS to global structure
    elseif detectionNullCondition && authorisedDetection 
        % DETECTION HAS BEEN LOST
        [METAObjUpdate,EVENT] = OMAS_eventHandler(SIM,METAObjUpdate,METAObjEvaluation,eventType.null_detection);
        metaEVENTS = vertcat(metaEVENTS,EVENT);                                                    % Add new META EVENTS to global structure
    end
    
    %% ASSESS WARNING CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    objectStatus = METAObjUpdate.objectStatus(index,eventType.warning);                            % Has a warning already been issued
    warningCondition = SIM.warningDistance >= euclideanSeperation && objectStatus == 0;            % The warning condition as a function of the object size
    warningNullCondition = SIM.warningDistance < euclideanSeperation && objectStatus == 1;         % The warning loss condition
    if warningCondition && isCollidable                                                                     
        % GENERATE THE WARNING EVENT
        [METAObjUpdate,EVENT] = OMAS_eventHandler(SIM,METAObjUpdate,METAObjEvaluation,eventType.warning);
        metaEVENTS = vertcat(metaEVENTS,EVENT);                                                    % Add new META EVENTS to global structure
    elseif warningNullCondition && isCollidable
        % GENERATE THE WARNING NULLIFICATION EVENT
        [METAObjUpdate,EVENT] = OMAS_eventHandler(SIM,METAObjUpdate,METAObjEvaluation,eventType.null_warning);
        metaEVENTS = vertcat(metaEVENTS,EVENT);                                                    % Add new META EVENTS to global structure
    end
    
    %% ASSESS COLLISION CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    objectStatus = METAObjUpdate.objectStatus(index,eventType.collision);                          % Has a collision already been issued
    collisionCondition = geometricLimits >= euclideanSeperation && objectStatus == 0;              % The collision condition as a function of the object size
    collisionNullCondition = geometricLimits < euclideanSeperation && objectStatus == 1;           % The conditions for breaking collision
    if collisionCondition && isCollidable
        % GENERATE THE COLLISION EVENT
        [METAObjUpdate,EVENT] = OMAS_eventHandler(SIM,METAObjUpdate,METAObjEvaluation,eventType.collision);
        metaEVENTS = vertcat(metaEVENTS,EVENT);                                                    % Add new META EVENTS to global structure
    elseif collisionNullCondition && isCollidable
        % GENERATE THE WARNING NULLIFICATION EVENT
        [METAObjUpdate,EVENT] = OMAS_eventHandler(SIM,METAObjUpdate,METAObjEvaluation,eventType.null_collision);
        metaEVENTS = vertcat(metaEVENTS,EVENT);                                                    % Add new META EVENTS to global structure
    end
    
    %% IF COLLIDABLE, NO WAYPOINT CHECK NECESSARY %%%%%%%%%%%%%%%%%%%%%%%%%
    if isCollidable
        SIM.OBJECTS(SIM.globalIDvector == METAObjUpdate.objectID) = METAObjUpdate;                          % Assign the updated SIM.OBJECT
        continue 
    end
    
    %% ASSESS WAYPOINT CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    objectStatus = METAObjUpdate.objectStatus(index,eventType.waypoint);                    
    waypointGetCondition = geometricLimits >= euclideanSeperation && objectStatus == 0;            % The collision condition as a function of the object size
    waypointNullCondition = geometricLimits < euclideanSeperation && objectStatus == 1;            % The conditions for breaking collision
    if waypointGetCondition && authorisedDetection
        % GENERATE THE WAYPOINT ACHIEVED EVENT
        [METAObjUpdate,EVENT] = OMAS_eventHandler(SIM,METAObjUpdate,METAObjEvaluation,eventType.waypoint);
        metaEVENTS = vertcat(metaEVENTS,EVENT);                                                    % Add new META EVENTS to global structure        
    elseif waypointNullCondition && authorisedDetection
        % GENERATE THE WAYPOINT NULLIFICATION EVENT
        % [METAObjUpdate,EVENT] = OMAS_eventHandler(SIM,METAObjUpdate,METAObjEvaluation,eventType.null_waypoint);
        % metaEVENTS = vertcat(metaEVENTS,EVENT);                                                  % Add new META EVENTS to global structure
    end
    
    %% UPDATE THE SIMULATION META COPY OF OBJECT %%%%%%%%%%%%%%%%%%%%%%%%%%  
    SIM.OBJECTS(SIM.globalIDvector == METAObjUpdate.objectID) = METAObjUpdate;                     % Update SIM structure for EVENTS and condition comparison
end
end
% UPDATE THE AGENT OBJECTS
function [referenceObject,objectEVENTS] = updateObjects(SIM,objectIndex,referenceObject)
% This function updates a referenceObject class against the rest of the
% objectIndex, independantly of the SIM.OBJECTS META data.
% INPUTS:
% SIM             - Current META structure
% objectIndex     - The current object list
% referenceObject - The agent object being updated

% OUTPUTS:
% referenceObject - The updated agent class

objectEVENTS = []; % Container for object based events

% GET THE OBJECTS EQUIVALENT SIM OBJECT
SIMreference = SIM.OBJECTS(SIM.globalIDvector == referenceObject.objectID);

if SIM.verboseMode
    fprintf('[UPDATE]\tCycling %s object (%s).\n',SIMreference.type,referenceObject.name);
end

% SWITCH BEHAVIOUR BASED ON OBJECT TYPE
switch SIMreference.type
    case OMAS_objectType.agent
        % AGENT - SIMULATION/ENVIROMENTAL FEEDBACK REQUIRED %%%%%%%%%%%%%%
        if SIM.verboseMode
            
        end
        observationPacket = []; % Container for agent detected objects
        
        % MOVE THROUGH THE COMPLETE OBJECT SET
        for METAindex = 1:SIM.totalObjects
            % GET THE SIM OBJECT OF THE EVALUATION OBJECT
            SIMevaluation = SIM.OBJECTS(METAindex);
            % GET DETECTION STATUS OF OBJECT
            objectStatus = SIMreference.objectStatus(METAindex,eventType.detection);   % Get the current status of this object
            
            % SKIP AGENT IF SELF COMPARISON OR CURRENT OBJECT UNDETECTED
            if SIMreference.objectID == SIMevaluation.objectID || 1 ~= objectStatus
                continue
            end
                            
            %% ELSE AGENT IS DETECTED (GENERATE A LOCALISED PROXIMITY DESCIPTION)
            % The second agent is within global seperation conditions. A
            % communication packet is therefore assembled containing the
            % data on the second object. rotated into the agent-local
            % coordinate frame.

            % GET THE RELATIVE GLOBAL POSITIONS OF THE AGENT AND OBJECT (x_rel = x_B - x_A)
            relativePosition = SIMreference.relativePositions(METAindex,:)';
            relativeVelocity = SIMevaluation.globalState(4:6) - SIMreference.globalState(4:6);
            % GET THE GLOBAL-(AGENT)LOCAL ROTATION MATRIX
            agent_R_GB = SIMreference.R_GB;
            % ROTATE THE GLOBAL STATE OF THE OBJECT INTO THE AGENT FRAME &
            % DISTORT THE OBSERVED POSITION (INTRODUCE PROBABILITY)
            % CARTESIAN REPRESENTATION
            observedPosition = (agent_R_GB*relativePosition); % Cartesian position
            observedVelocity = (agent_R_GB*relativeVelocity); % Cartesion velocity            
            % SPHERICAL REPRESENTATION
            observedRange     = norm(observedPosition);
            observedElevation = asin(observedPosition(3)/observedRange);
            observedAzimuth   = atan2(observedPosition(2),observedPosition(1));
            % PRESENT SIZE PROPERTIES
            observedSize      = objectIndex{METAindex}.VIRTUAL.size;  % Size
            observedAngularWidth = 2*asin(observedSize/(observedRange + observedSize));               
            
            % OBJECT TIME TO COLLISION APPROXIMATION
            tau = -(dot(observedPosition,observedVelocity)/dot(observedVelocity,observedVelocity)); % Projection based approximation
            if isnan(tau)
                tau = inf;
            end            
            
            % OBJECT PRIORITY (IF WAYPOINT)
            if SIMevaluation.type == OMAS_objectType.waypoint
                association = objectIndex{METAindex}.getAgentAssociation(referenceObject);
                observedPriority = association.priority;
            else
                observedPriority = NaN;
            end
            
            %% PREPARE DETECTION PACKET BASED ON LOCAL DATA (FROM GLOBAL)
            detectionObject = struct('objectID',SIMevaluation.objectID,...     % The object ID (observed)
                                     'name',SIMevaluation.name,...             % The object name tag (observed)
                                     'type',SIMevaluation.type,...             % The objects sim-type enum
                                     'timeToCollision',tau,...                 % The objects geometric time to collision (velocity vector comparision)
                                     'priority',observedPriority,...           % The objects global priority for that agent
                                     'size',observedSize,...                   % The objects true size
                                     'angularWidth',observedAngularWidth,...   % The apparent angular width at that range
                                     'position',observedPosition,...           % The apparent position in the relative frame
                                     'velocity',observedVelocity,...           % The apparent velocity in the relative frame
                                     'range',observedRange,...                 % The apparent range
                                     'inclinationAngle',observedElevation,...  % The apparent inclination angle
                                     'azimuthAngle',observedAzimuth,...        % The apparent Azimuth angle
                                     'colour',SIMevaluation.colour);           % Finally the simulation's colourID          
            observationPacket = vertcat(observationPacket,detectionObject);    % Append object to packet to agent
        end
        
        % COMPUTE THE LOCAL AGENT CYCLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Given the detection object defined for this agent, compute the
        % agents cycle with this information.
%         try
            % SEND OBJECT OBSERVERATION PACKET TO AGENT
            referenceObject = referenceObject.processTimeCycle(SIM.TIME,observationPacket);                      % Hand the current sim TIME object and complete observation packet
%         catch objectProcessError                                                                                 % to the reference object/agent
%             warning('[ERROR]\tProcess time cycle for agent %s',referenceObject.name);
%             error(objectProcessError.message);
%         end
    case OMAS_objectType.obstacle
        % PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        referenceObject = referenceObject.processTimeCycle(SIM.TIME);
    case OMAS_objectType.waypoint
        % WAYPOINTS ARE CONSIDERED PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%
        referenceObject = referenceObject.processTimeCycle(SIM.TIME);
    otherwise
        % CONSIDERED PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%%%%%%%%%%%%%%%
        referenceObject = referenceObject.processTimeCycle(SIM.TIME);
end
end

%% SIMULATION OUTPUT OPERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET INITIAL OUTPUT DATA STRUCTURE
function [DATA] = getOutputStructure(SIM)
% This function assembles the initial output DATA structure from the
% simulation input variables and it back as a simulation output container.
% INPUT:
% SIM  - A local reference to the META structure
% OUPUT:
% DATA - The output data structure

% CALCULATE THE INITIAL DATA PARAMETERS
globalStateNum = size(SIM.OBJECTS(1).globalState,1);
systemStates = globalStateNum*SIM.totalObjects;                            % The total number of system states
% DEFINE THE GLOBAL STATE INDICES
indexSet = zeros(2,length(SIM.OBJECTS));                                   % Indices are defined [start;end]*n
indexSet(1,:) = (0:globalStateNum:systemStates-globalStateNum) + 1;        % Declare the state indices
indexSet(2,:) = globalStateNum:globalStateNum:systemStates;

% GENERATE THE OUTPUT STRUCTURE FROM THE SIMULATION INPUT PARAMETERS
DATA = struct('outputPath',strcat(SIM.outputPath,'DATA.mat'),...
               'outputDir',SIM.outputPath,...
            'totalObjects',SIM.totalObjects,...
             'totalAgents',SIM.totalAgents,...
          'totalObstacles',SIM.totalObstacles,...
          'totalWaypoints',SIM.totalWaypoints,...
              'timeVector',SIM.TIME.timeVector,...
                      'dt',SIM.TIME.dt,...
              'stateIndex',indexSet,...
      'globalTrajectories',zeros(systemStates,SIM.TIME.numSteps));         % Prepare the output container
end
% SAVE DATA TO FILES
function sendOutputToFiles(META,EVENTS,DATA,objectIndex)
% This function is designed to handle the output data from the simulation
% and export the variables to background files
% INPUTS:
% META        - Local copy of the META variable
% EVENTS      - The cell array of event history objects
% DATA        - The output data structure
% objectIndex - The object (non-META) object class cell array

% DISPLAY ALL VARIABLES
whos;
% SAVE META DATA
save(strcat(META.outputPath,'META.mat'),'META');
save(strcat(META.outputPath,'OBJECTS.mat'),'objectIndex')
% SAVE EVENT HISTORY
save(strcat(META.outputPath,'EVENTS.mat'),'EVENTS');
% SAVE OUTPUT DATA
save(strcat(META.outputPath,'DATA.mat'),'DATA');
% DISPLAY NOTIFICATION
fprintf('[%s]\tData objects outputted to file:\n',META.phase);
fprintf('[%s]\tDirectory: %s\n',META.phase,META.outputPath);
end
