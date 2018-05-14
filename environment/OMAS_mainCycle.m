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

%% PRE-SIMULATION DECLARATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
META.phase = 'META'; % Declare the simuation active phase
EVENTS = [];         % Reset the EVENT log container

% PREPARE THE OUTPUT DATA CONTAINER
[DATA] = getOutputStructure(META); 

%%%%%%%%%%%%%%%%%%%%%% BEGIN TIME-STEP INTERATIONS %%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[%s]\n[%s]\tLAUNCHING SIMULATION...\n[%s]\n',META.phase,META.phase,META.phase);
% EXECUTE OMAS MAIN CYCLE
[META,objectIndex,DATA,EVENTS] = processTimeSeries(META,objectIndex,DATA,EVENTS);
% DATA.computationTime = toc; % Get total elapsed simulation time
% fprintf('\n[%s]\tSIMULATION COMPLETE (time elapsed: %ss)\n',META.phase,num2str(DATA.computationTime));
fprintf('[%s]\n[%s]\t...SIMULATION COMPLETE\n[%s]\n',META.phase,META.phase,META.phase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FINALISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
META.phase = 'OUTPUT';                                   % Change to OUTPUT phase
% clearvars -except META EVENTS DATA objectIndex         % Clear loose simulation variables
fprintf('[%s]\tDumping initial simulation variables to output directory\n[%s]\n',META.phase,META.phase);
% sendOutputToFiles(META,EVENTS,DATA,objectIndex);
fprintf('[%s]\tClosing simulation.\n',META.phase);
SIM = META;
end

%%%%%%%%%%%%%%%%%%%%%%%% MAIN CYCLE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% for step = 1:META.TIME.numSteps
step = 1; 
while step <= META.TIME.numSteps
    
    % 0. /////////////////// CHECK AGENT IDLE STATUS //////////////////////
    % If there are agents in the simulation, then break if all are idle,
    % else allow run to full duration
    if sum([META.OBJECTS([META.OBJECTS.type] == OMAS_objectType.agent).idleStatus]) == META.totalAgents && META.totalAgents > 0
        fprintf('[%s]\t... All agents are idle, abort.\n',META.phase);
        break
    end
    
    % 1. ////////////// UPDATE TIMESTEP PARAMETERS (@t=k) ////////////////
    META.TIME.currentTime = META.TIME.timeVector(step);                    % The current sim-time
    META.TIME.currentStep = step;                                          % The current sim-step
    if META.verbosity >= 1
        fprintf('[%s]\tStep: %.0f\tTime: %.2fs\n',META.phase,META.TIME.currentStep,META.TIME.currentTime);
    end
    
    % 2. ///////////////// RECORD GLOBAL STATE (@t=k) ///////////////////////
    for ID1 = 1:META.totalObjects
        % Collect the META.OBJECTS.state data (as the simulations understanding of the global states) 
        % and save to the output DATA.globalTrajectories set, must be done synchronously).
        DATA.globalTrajectories(DATA.stateIndex(1,ID1):DATA.stateIndex(2,ID1),META.TIME.currentStep) = META.OBJECTS(ID1).globalState;
    end

    % 3. ////////////// UPDATE THE GLOBAL STATES (@t=k) ////////////////////
    for ID1 = 1:META.totalObjects                                          % For each META object element
        META.OBJECTS(ID1) = updateGlobalStates(META,objectIndex{ID1}.objectID,objectIndex{ID1}.VIRTUAL);
    end
    
    % 4. ////////// UPDATE THE VISUAL REPRESENTATIONS (@t=k) //////////////
    for ID1 = 1:META.totalObjects
        if ~isstruct(META.OBJECTS(ID1).patch)
            continue
        end
        % MODIFY THE ASSOCIATED GEOMETRY
        META.OBJECTS(ID1) = updateVisuals(META.OBJECTS(ID1));
    end
    
    % 5. //////// UPDATE SIMULATION/ENVIRONMENTAL META DATA (@t=k) ////////
    objectSnapshot = objectIndex;                                          % Make a temporary record of the object set
    for ID1 = 1:META.totalObjects                                          % For each META object element
        % UPDATE EACH META.OBJECT AGAINST THE CURRENT META AND OBJECT 
        % INDEX SNAPSHOT
        [META.OBJECTS(ID1),metaEVENTS] = updateSimulationMeta(META,objectSnapshot,META.OBJECTS(ID1));  % Update META snapshot with equivilant objectIndex.state
        % LOG THE META EVENTS
        if ~isempty(metaEVENTS)                                            % If META events occurred this timestep
            EVENTS = vertcat(EVENTS,metaEVENTS);                           % Append to global EVENTS
        end
    end
    
    % 6. //// COMPUTE AGENT CYCLES (UPDATE GLOBAL REPRESENTATION (@t=k) ///
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
    % 7. /// THE 'OBJECT.VIRTUAL' PROPERTIES IS NOW UPDATED FOR (t=k+1) ///
    step = step + 1;
end
% CREATE TERMINAL VALUES FOR CLARITY
META.TIME.endStep = META.TIME.currentStep;
META.TIME.endTime = META.TIME.currentTime; 
end

%%%%%%%%%%%%%%%%%%%%%%%%% UPDATE PROCEDURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE THE GLOBAL STATE VECTORS FROM ENTITY.VIRTUAL
function [METAObjUpdate] = updateGlobalStates(SIM,objectID,objectVIRTUALparameters)
% This function reallocates the global properties of each entity to the 
% META.OBJECT series to allow faster reference and increased independance
% of the main cycle from the object cycles.

% IDENTIFY THE ASSOCIATED META OBJECT 
METAObjUpdate  = SIM.OBJECTS(SIM.globalIDvector == objectID);              % Get the current META object associated with 'entity'

% EXTRACT THE UPDATED GLOBAL PARAMETERS FOR THE ENTITY
globalToRotatedBodyQuaternion_k = objectVIRTUALparameters.quaternion;                      % Get the global triad to rotated body rotation
globalVelocity_k = objectVIRTUALparameters.globalVelocity;

% CONFIRM DIMENSIONS
assert(size(globalToRotatedBodyQuaternion_k,1) == 4,'Object quaternion update must be given as a column vector [4x1].');
assert(size(globalVelocity_k,1) == 3,'Object velocity update must be given as a column vector [3x1].');

% //////////////////////// UPDATE META.OBJECT /////////////////////////////

% UPDATE THE GLOBAL TO BODY & BODY TO GLOBAL ROTATION MATRIX FROM THE NEW 
[METAObjUpdate.R_BG,METAObjUpdate.R_GB] = OMAS_axisTools.quaternionToRotationMatrix(globalToRotatedBodyQuaternion_k);
% UPDATE THE GLOBAL POSITION
globalPosition_k = METAObjUpdate.globalState(1:3) + globalVelocity_k*SIM.TIME.dt; % Calculate the new global position
% REBUILD GLOBAL STATES (ENTITY & META)
METAObjUpdate.globalState = [globalPosition_k;globalVelocity_k;globalToRotatedBodyQuaternion_k]; 

% DETERMINE IF ENTITY HAS INDICATED THAT TASK IS COMPLETE
METAObjUpdate.idleStatus = objectVIRTUALparameters.idleStatus;
end
% UPDATE THE OBJECT META PROPERTIES
function [METAObjUpdate,metaEVENTS] = updateSimulationMeta(SIM,objectProofSet,METAObjUpdate)
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

% REFERENCE OBJECT CONSTANTS
metaEVENTS = [];                                                           % Reset meta Events container
detectionRange = SIM.visabilityModifier*(METAObjUpdate.detectionRange);    % Get the detection range from newMETAObj

% ITERATE THROUGH THE SIMULATION OBJECT SET
for index = 1:SIM.totalObjects                                             % Loop each agent record in the META structure
    % UPDATE GLOBAL PROPERTIES (META DATA)
    % GET THE OBJECT TO EVALUATE AGAINST "METAObjUpdate" 
    METAObjEvaluation = SIM.OBJECTS(index);
    % NEGLECT SELF EVALUATION
    if METAObjEvaluation.objectID == METAObjUpdate.objectID                % Where not own ID1
       continue 
    end

    % THE NEW GLOBAL STATE (METAObjUpdate) IS NOW RE-EVALUATED AGAINST THE OTHER SIM.OBJECT DATA
    % Define the new relative global state (x_rel = x_b - x_a)
    relativePosition = METAObjEvaluation.globalState(1:3) - METAObjUpdate.globalState(1:3);          % Relative global position
    METAObjUpdate.relativePositions(index,:) = relativePosition;                                     % Allocate new obj relative position vector
    euclideanSeparation = norm(relativePosition);                                                    % Scalar seperation value
    METAObjUpdate.euclideanSeparation(index,1) = euclideanSeparation;                                % Get the seperation for validation
    geometricLimits = METAObjEvaluation.radius + METAObjUpdate.radius;                               % Evaluate the combined agent radii
    
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
    detectionCondition = detectionRange >= euclideanSeparation && objectStatus == 0;               % Detection condition (all objects)
    detectionNullCondition = detectionRange < euclideanSeparation && objectStatus == 1;            % Detection loss condition (all objects)
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
    warningCondition = SIM.warningDistance >= euclideanSeparation && objectStatus == 0;            % The warning condition as a function of the object size
    warningNullCondition = SIM.warningDistance < euclideanSeparation && objectStatus == 1;         % The warning loss condition
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
    collisionCondition = euclideanSeparation < (geometricLimits - SIM.conditionTolerance) && objectStatus == 0; % The collision condition as a function of the object size
    collisionNullCondition = geometricLimits < euclideanSeparation && objectStatus == 1;           % The conditions for breaking collision
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
    waypointGetCondition = geometricLimits >= euclideanSeparation && objectStatus == 0;            % The collision condition as a function of the object size
    waypointNullCondition = geometricLimits < euclideanSeparation && objectStatus == 1;            % The conditions for breaking collision
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
% UPDATE THE OBJECT PROPERTIES
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

if SIM.verbosity >= 2
    fprintf('[UPDATE]\tCycling %s\t(%s).\n',char(SIMreference.type),referenceObject.name);
end

% SWITCH BEHAVIOUR BASED ON OBJECT TYPE
switch SIMreference.type
    case OMAS_objectType.agent
        % AGENT - SIMULATION/ENVIROMENTAL FEEDBACK REQUIRED %%%%%%%%%%%%%%
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
            % ROTATE THE GLOBAL STATE OF THE OBJECT INTO THE AGENT FRAME &
            % DISTORT THE OBSERVED POSITION (INTRODUCE PROBABILITY)
            % CARTESIAN REPRESENTATION
            observedPosition = SIMreference.R_GB*relativePosition;         % Cartesian position
            observedVelocity = SIMreference.R_GB*relativeVelocity;         % Cartesion velocity            
            % SPHERICAL REPRESENTATION
            observedRange     = norm(observedPosition);
            observedElevation = asin(observedPosition(3)/observedRange);
            observedAzimuth   = atan2(observedPosition(2),observedPosition(1));
            % PRESENT SIZE PROPERTIES
            observedSize      = objectIndex{METAindex}.VIRTUAL.radius;     
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
                                     'radius',observedSize,...                 % The objects true size
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
%         catch objectProcessError 
%             warning(objectProcessError.message);
%             error('OpenMas: Cycle failed, there was a problem with the object (%s) file, please check: %s',...
%                   referenceObject.name,objectProcessError.stack(1).name);    
%         end
    case OMAS_objectType.obstacle
        % PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        referenceObject = referenceObject.processTimeCycle(SIM.TIME);
    case OMAS_objectType.waypoint
        % WAYPOINTS ARE CONSIDERED PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%
        referenceObject = referenceObject.processTimeCycle(SIM.TIME);
    otherwise
        % CONSIDERED PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%%%%%%%%%%%%%%%
        % DO NOT UPDATE
end
end
% UPDATE VISUAL REPRESENTATION OF THE OBJECTS (IF EXISTS)
function [METAobject] = updateVisuals(METAobject)
% This utility preforms the visualisation update for OpenMAS as the time
% progression continues.

% OBJECT PROPERTIES
% position = METAobject.globalState(1:3);
% velocity = METAobject.globalState(4:6);
% quaternion = METAobject.globalState(7:10);
% eulerRotations = OMAS_axisTools.quaternionToEuler(quaternion);

%verts = METAobject.patch.vertices*METAobject.R_BG + METAobject.globalState(1:3)';

% % CONVERT TO VR AXIS CONVENTION
% positionVR = position;
% positionVR(2) = -position(3);
% positionVR(3) = position(2);        % The objects global VR position
% velocityVR = velocity;
% velocityVR(2) = -velocity(3);
% velocityVR(3) = velocity(2);        % The objects global VR velocity
% rotationVR = zeros(3,1);
% rotationVR(1) = -eulerRotations(2);
% rotationVR(2) = eulerRotations(3);
% rotationVR(3) = -eulerRotations(1); % The VR rotations from the body rotations
% % GET THE DCM DESCRIBING THE ATTITUDE UPDATE
% DCM_update = simulation_axisTools.eulerToRotationMatrix(rotationVR);

% FROM THE DCM MATRIX, GET THE ROTATED VR-MODEL ROTATION

% display('UPDATING VR OBJECT');

% THE OBJECTS VRML MUST BE RELATED TO THE META.OBJECT.objectID
end

%%%%%%%%%%%%%%%%%%%% SIMULATION OUTPUT OPERATIONS %%%%%%%%%%%%%%%%%%%%%%%%%
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
DATA = struct('outputPath',[SIM.outputPath,'DATA.mat'],...
               'outputDir',SIM.outputPath,...
            'totalObjects',SIM.totalObjects,...
             'totalAgents',SIM.totalAgents,...
          'totalObstacles',SIM.totalObstacles,...
          'totalWaypoints',SIM.totalWaypoints,...
              'timeVector',SIM.TIME.timeVector,...
                      'dt',SIM.TIME.dt,...
              'stateIndex',indexSet,...
      'globalTrajectories',NaN(systemStates,SIM.TIME.numSteps));         % Prepare the output container
%       'globalTrajectories',zeros(systemStates,SIM.TIME.numSteps));         % Prepare the output container
end