%% THE OPENMAS CYCLE PROCESSOR (OMAS_process.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the progressive cycles of the simulation, for the
% environment, obstacles, agents and objects withing the simulation space.

% Author: James A. Douthwaite 06/10/2016

% //////////////////////// MAIN WRAPPER SCRIPT ////////////////////////////
function [DATA,SIM,EVENTS,objectIndex]  = OMAS_process(META,objectIndex)
% INPUTS:
% objectIndex - The cell array of object classes
% OUTPUTS:
% DATA        - The output DATA structure with timeseries data
% SIM         - The terminal META structure copy

%% PRE-SIMULATION DECLARATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
META.phase = 'META'; % Declare the simuation active phase
EVENTS = [];         % Reset the EVENT log container

% PREPARE THE OUTPUT DATA CONTAINER
[DATA] = GetOutputStructure(META); 

%%%%%%%%%%%%%%%%%%%%%% BEGIN TIME-STEP INTERATIONS %%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[%s]\n[%s]\tLAUNCHING SIMULATION...\n[%s]\n',META.phase,META.phase,META.phase);
% EXECUTE OMAS MAIN CYCLE
[META,objectIndex,DATA,EVENTS] = ProcessGlobalTimeSeries(META,objectIndex,DATA,EVENTS);
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

% //////////////////////// MAIN CYCLE FUNCTIONS ///////////////////////////
% PROCESS TIME VECTOR
function [META,objectIndex,DATA,EVENTS] = ProcessGlobalTimeSeries(META,objectIndex,DATA,EVENTS)
% This function computes the simulation loop cycle across the given
% descrete META.timevector:
% INPUTS:
% META - The initialised META data-structure
% objectIndex - The cell vector of objects (and/or agents)
% DATA        - The initialised output DATA structure
% EVENTS      - The initialised EVENTS structure list

% PROCESS TIME CYCLE
step = META.TIME.currentStep; 
exitFlag = 0; 
while step <= META.TIME.numSteps
    
    % 0. /////////////////// CHECK AGENT IDLE STATUS //////////////////////
    % If there are agents in the simulation, then break if all are idle,
    % else allow run to full duration
    idleCondition = sum([META.OBJECTS([META.OBJECTS.type] == OMAS_objectType.agent).idleStatus]) == META.totalAgents && META.totalAgents > 0;
    if idleCondition == 1 && ~exitFlag
        exitFlag = 1; 
        referenceTime = META.TIME.currentTime;
        fprintf('[%s]\t... All agents are idle, aborting in t=%0.2fs.\n',META.phase,META.TIME.idleTimeOut);
    elseif exitFlag == 1
        timeOut = META.TIME.idleTimeOut - (META.TIME.currentTime - referenceTime);
        if timeOut <= 0
           fprintf('[%s]\t... Aborting.\n',META.phase); 
           break 
        end        
        fprintf('[%s]\t... All agents are idle, aborting in t=%0.2fs.\n',META.phase,timeOut);
    end 
    % /////////////////////////////////////////////////////////////////////
    
    % 1. ////////////// UPDATE TIMESTEP PARAMETERS (@t=k) /////////////////
    META.TIME.currentTime = META.TIME.timeVector(step);                    % The current sim-time
    META.TIME.currentStep = step;                                          % The current sim-step
    if META.verbosity >= 1
        fprintf('[%s]\tStep: %.0f\tTime: %.2fs\n',META.phase,META.TIME.currentStep,META.TIME.currentTime);
    end  
    % /////////////////////////////////////////////////////////////////////
    
    % 2. ///////////////// RECORD GLOBAL STATE (@t=k) /////////////////////
    for ID1 = 1:META.totalObjects
        % Collect the META.OBJECTS.state data (as the simulations understanding of the global states) 
        % and save to the output DATA.globalTrajectories set, must be done synchronously).
        DATA.globalTrajectories(DATA.stateIndex(1,ID1):DATA.stateIndex(2,ID1),META.TIME.currentStep) = META.OBJECTS(ID1).globalState;
    end
    % /////////////////////////////////////////////////////////////////////
    
    % 3. ////////////// UPDATE THE GLOBAL STATES (@t=k) ////////////////////
    for ID1 = 1:META.totalObjects                                          % For each META object element
        % Update META.OBJECT the structure
%         META.OBJECTS(ID1) = OMAS_updateGlobalStates_mex(...
%             META,...
%             objectIndex{ID1}.objectID,...
%             objectIndex{ID1}.GetGLOBAL('velocity'),...
%             objectIndex{ID1}.GetGLOBAL('quaternion'),...
%             objectIndex{ID1}.GetGLOBAL('idleStatus'));   
        META.OBJECTS(ID1) = OMAS_updateGlobalStates(...
            META,...
            objectIndex{ID1}.objectID,...
            objectIndex{ID1}.GetGLOBAL('velocity'),...
            objectIndex{ID1}.GetGLOBAL('quaternion'),...
            objectIndex{ID1}.GetGLOBAL('idleStatus')); 
    end
    % /////////////////////////////////////////////////////////////////////
    
%     % 4. ////////// UPDATE THE VISUAL REPRESENTATIONS (@t=k) //////////////
%     for ID1 = 1:META.totalObjects
%         if ~isstruct(META.OBJECTS(ID1).geometry)
%             continue
%         end
%         % MODIFY THE ASSOCIATED GEOMETRY
%         META.OBJECTS(ID1) = updateVisuals(META.OBJECTS(ID1));
%     end
    
    % 5. //////// UPDATE SIMULATION/ENVIRONMENTAL META DATA (@t=k) ////////
    [META,metaEVENTS] = UpdateSimulationMeta(META,objectIndex);            % Update META snapshot with equivilant objectIndex.state
    % LOG THE META EVENTS
    if ~isempty(metaEVENTS)                                                % If META events occurred this timestep
        EVENTS = vertcat(EVENTS,metaEVENTS);                               % Append to global EVENTS
    end
    % /////////////////////////////////////////////////////////////////////
    
    % 6. //// COMPUTE AGENT CYCLES (UPDATE GLOBAL REPRESENTATION (@t=k) ///
    if META.threadPool ~= 0
        objectSnapshot = objectIndex;                                      % Make a temporary record of the object set
        parfor (ID1 = 1:META.totalObjects)
            % MOVE THROUGH OBJECT INDEX AND UPDATE EACH AGENT
            [objectIndex{ID1},objectEVENTS] = UpdateObjects(META,objectSnapshot,objectIndex{ID1});            
            % LOG THE OBJECT EVENTS
            if ~isempty(objectEVENTS)
                EVENTS = vertcat(EVENTS,objectEVENTS);
            end
        end
    else
        for ID1 = 1:META.totalObjects
            % MOVE THROUGH OBJECT INDEX AND UPDATE EACH AGENT
            [objectIndex{ID1},objectEVENTS] = UpdateObjects(META,objectIndex,objectIndex{ID1}); % Update objectIndex snapshot with new META data
            % LOG THE OBJECT EVENTS
            if ~isempty(objectEVENTS)                                      % If objectEVENTS occur in this timestep
                EVENTS = vertcat(EVENTS,objectEVENTS);                     % Append to global EVENTS
            end
        end
    end
    % /////////////////////////////////////////////////////////////////////
    
    % 7. /// THE 'OBJECT.VIRTUAL' PROPERTIES IS NOW UPDATED FOR (t=k+1) ///
    step = step + 1;
end
% CREATE TERMINAL VALUES, FOR CLARITY
META.TIME.endStep = META.TIME.currentStep;
META.TIME.endTime = META.TIME.currentTime; 
end
% UPDATE THE OBJECT META PROPERTIES
function [SIM,metaEVENTS]	= UpdateSimulationMeta(SIM,objectIndex)
% In this function we wish to update the global separations and event
% conditions based on the new 'globalState' properties.

% INPUTS:
% SIM         - A local copy of the global META structure
% objectIndex - The cell array of object classes.
% OUTPUTS:
% SIM         - The updated META structure
% metaEVENTS  - A vector of new event objects.

% We want to move through the object set an update both the object being
% updated, but also the object it is being updated against.

% //// ASSESS THE COLLISIONS AND RELATIVE POSITIONS (SYMMETRICAL CHECK) ///
collisionLogicals = zeros(SIM.totalObjects);
warningLogicals   = zeros(SIM.totalObjects);
for entityA = 1:SIM.totalObjects                                           % Object A's position in the META.OBJECTS
    
    % NOTE:
    % - We only need to consider the objects whos separations have yet to be
    %   evaluated. If we update 1&2 and 2&1 simultaneously, then we must only
    %   examine the unique ID permutations to assess their conditions.
    
    for entityB = (entityA+1):(SIM.totalObjects)                           % Object B's position in the META.OBJECTS
        % THE OBJECTS ASSOCIATED objectIndex
        object_A = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entityA).objectID};
        object_B = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entityB).objectID};  % Determine the associated geometries
        
        % UPDATE RELATIVE POSITIONS  (centroid separation based)
        SIM.OBJECTS(entityA).relativePositions(entityB,:) = (SIM.OBJECTS(entityB).globalState(1:3,1) - SIM.OBJECTS(entityA).globalState(1:3,1))';   % AB vector
        SIM.OBJECTS(entityB).relativePositions(entityA,:) = (SIM.OBJECTS(entityA).globalState(1:3,1) - SIM.OBJECTS(entityB).globalState(1:3,1))';   % BA vector
        
        % GET THE WARNING CONDITIONS (centroid separation based)
        separationDistance = norm(SIM.OBJECTS(entityA).relativePositions(entityB,:)) - (SIM.OBJECTS(entityA).radius + SIM.OBJECTS(entityB).radius);
        warningCondition = SIM.warningDistance >= separationDistance;
        warningLogicals(entityA,entityB) = warningCondition;
        warningLogicals(entityB,entityA) = warningCondition;               % Assign the warning condition
        
        % EVALUATE COLLISIONS BETWEEN THE OBJECTS
        collisionCondition = OMAS_collisionDetection(...
            SIM.OBJECTS(entityA),object_A.GEOMETRY,...
            SIM.OBJECTS(entityB),object_B.GEOMETRY,...
            SIM.conditionTolerance);
        collisionLogicals(entityA,entityB) = collisionCondition;
        collisionLogicals(entityB,entityA) = collisionCondition;           % Assign the collision condition     
        
        % ERROR CHECKING (Should never happen)
        if separationDistance > 0 && collisionCondition
            error('[ERROR] A collision occurred at distance %f without violating the minimum separation of %f',...
                    norm(SIM.OBJECTS(entityA).relativePositions(entityB,:)),(SIM.OBJECTS(entityA).radius + SIM.OBJECTS(entityB).radius));
        end
    end
end

% //////////// ASSESS DETECTION CONDITIONS (ASYMMETRICAL CHECK) ///////////
% DEFAULT DETECTION CONDITION
detectionLogicals = zeros(SIM.totalObjects);
for entityA = 1:SIM.totalObjects
    % Check objects have the capacity to observe other objects
    if SIM.OBJECTS(entityA).type ~= OMAS_objectType.agent
        continue
    end
    
    % Check all agent observations
    for entityB = 1:SIM.totalObjects                                       % Object B's position in the META.OBJECTS
        % SKIP SELF-REFERENCE    
        if SIM.OBJECTS(entityA).objectID == SIM.OBJECTS(entityB).objectID                                           
            continue                                                       % Only agents can generate notifications
        end 
        % NOTE:
        % - We also need to determine which agents are able to observe other
        %   agents in the field. However, this is NOT A SYMMETRIC operation
        %   as A may be able to detect B without B detecting A. 
        
        % Detection check #1 - Authorisation
        isWaypoint = SIM.OBJECTS(entityB).type == OMAS_objectType.waypoint;
        if isWaypoint
            isAuthorised = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entityB).objectID}.IDAssociationCheck(SIM.OBJECTS(entityA).objectID);
            if ~isAuthorised
                continue
            end
        end
        
        % Detection check #2 - Spherical 
        isDetected = OMAS_geometry.intersect_spheres(...
            SIM.OBJECTS(entityA).globalState(1:3),SIM.OBJECTS(entityA).detectionRadius + 0.5*SIM.conditionTolerance,...
            SIM.OBJECTS(entityB).globalState(1:3),SIM.OBJECTS(entityB).radius + 0.5*SIM.conditionTolerance);
        
        % Object is not detected
        if ~isDetected
            continue
        end
        
        % This check is sufficient for way-points
        if isWaypoint
            detectionLogicals(entityA,entityB) = 1;
            continue
        end
        
        % ///////////////// OBJECT MAY BE OBSERVED ////////////////////////
        % NOTE:
        % - If the distance(less the radius) is less than detection radius
        %   then its worth checking for more a more complex representation.
        % - We want to check if any vertices/edges can be observed. If any 
        %   edges are visible, then they are to be sent to the first object.
        
        % Get the geometry of the second object (GEOMETRY OF B)
        geometryB = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entityB).objectID}.GEOMETRY;        
         
        % Check if a better check is necessary
        hasGeometry = isstruct(geometryB) && size(geometryB.vertices,1) > 0;    
        if isDetected && ~hasGeometry 
            detectionLogicals(entityA,entityB) = 1;                        % Detection occurred
            break
        end
        isDetected = 0; % Reset, potential false positive when checking against polygon
         
        % Assess the polygon
        relativeR = SIM.OBJECTS(entityA).R'*SIM.OBJECTS(entityB).R;        % Rotate second geometry orientated relative to
        % Design the vertex set
        geometryB.vertices = geometryB.vertices*relativeR + SIM.OBJECTS(entityB).globalState(1:3)';  
        
        % Assess each face and vertex for detection
        for face = 1:size(geometryB.vertices,1)
            % GET THE MEMBER IDs OF THE FACE
            faceMembers  = geometryB.faces(face,:);
            faceVertices = geometryB.vertices(faceMembers,:);
            
            % CHECK THE GEOMETRY AGAINST THE SPHERICAL CONSTRAINT (with mex)
%             isDetected = CheckSphereTriangleIntersection_mex(...
%                 SIM.OBJECTS(entityA).globalState(1:3),...
%                 SIM.OBJECTS(entityA).detectionRadius,...
%                 faceVertices(1,:)',faceVertices(2,:)',faceVertices(3,:)');
            % CHECK THE GEOMETRY AGAINST THE SPHERICAL CONSTRAINT
            isDetected = CheckSphereTriangleIntersection(...
                SIM.OBJECTS(entityA).globalState(1:3),...
                SIM.OBJECTS(entityA).detectionRadius,...
                faceVertices(1,:)',faceVertices(2,:)',faceVertices(3,:)');
            
            % IF THE FACE IS DETECTED
            if isDetected
                detectionLogicals(entityA,entityB) = 1;
                break
            end
        end
        
        % NOTE:
        % - If any faces are within the constraint, then the detection
        %   logical is set true and then exited.
        % - Otherwise, if the none of the contraints are met, then it
        %   remains undetected.
    end
end

% NOTES:
% - The fields geometry based fields of SIM.OBJECTS are updated to this
%   time-step. The event conditions can now be re-evaluted.
% - CHECKS MUST BE ASYMMETRICAL - A may detect B without B detecting A

metaEVENTS = [];                                                           % Reset meta Events container
for entityA = 1:SIM.totalObjects
    % CHECK OBJECTS HAVE THE CAPACITY TO OBSERVE OTHER OBJECTS
    if SIM.OBJECTS(entityA).type ~= OMAS_objectType.agent
        continue
    end
    
    for entityB = 1:SIM.totalObjects
        % OMIT SELF-CHECK
        if SIM.OBJECTS(entityA).objectID == SIM.OBJECTS(entityB).objectID                                         
            continue                                                       % Only agents can generate notifications
        end 
                
        % //////////// ASSESS EVENT CONDITIONS FOR THE AGENTS /////////////
        % NOTES:
        % - The agent set is only capable of generating events associated
        %   with the detections, warnings and collisions.
        
        % ///////////////// ASSESS DETECTIONS CONDITIONS //////////////////
        % DETECTION EVENT LOGIC
        ABdetectionConstraint       =  detectionLogicals(entityA,entityB);
        novelDetectionConstraint    =  SIM.OBJECTS(entityA).objectStatus(entityB,eventType.detection) == 0;
        detectionEventCondition     =  ABdetectionConstraint &&  novelDetectionConstraint;   % Is within a range and a detection event has not been issued.
        detectionEventNullCondition = ~ABdetectionConstraint && ~novelDetectionConstraint;   % Is outside a range and a detection condition is still toggled.
        % EVALUATE OCCURANCE
        if detectionEventCondition 
            % UPDATE DETECTION STATUS
            SIM.OBJECTS(entityA).objectStatus(entityB,eventType.detection) = 1;              % Ammend status of the META object
            % GENERATE DETECTION EVENT
            [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.detection);
            metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
        elseif detectionEventNullCondition
            % AMEND THE META OBJECT
            SIM.OBJECTS(entityA).objectStatus(entityB,eventType.detection) = 0;              % Ammend status of the META object
            % GENERATE THE DETECTION NULLIFICATION EVENT
            [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.null_detection);
            metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
        end
                
        % ///////////////// ASSESS WAYPOINT CONDITIONS ////////////////////
        if SIM.OBJECTS(entityB).type == OMAS_objectType.waypoint
            % Check A is allowed to get B
            waypointObj = objectIndex{SIM.globalIDvector == SIM.OBJECTS(entityB).objectID};
            hasAssociation = waypointObj.IDAssociationCheck(SIM.OBJECTS(entityA).objectID);
            % WAYPOINT EVENT LOGIC
            ABwaypointCondition        =  collisionLogicals(entityA,entityB) && hasAssociation;
            novelWaypointGetConstraint =  SIM.OBJECTS(entityA).objectStatus(entityB,eventType.waypoint) == 0;
            waypointGetEventCondition  =  ABwaypointCondition &&  novelWaypointGetConstraint;    % Is within a range and a collision event has not been issued.
            waypointEventNullCondition = ~ABwaypointCondition && ~novelWaypointGetConstraint;    % Is outside a range and a waypoint-get condition is still toggled.
            % EVALUATE OCCURANCE
            if waypointGetEventCondition
                % AMEND THE META OBJECT
                SIM.OBJECTS(entityA).objectStatus(entityB,eventType.waypoint) = 1;
                % GENERATE THE WAYPOINT ACHIEVED EVENT
                [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.waypoint);
                metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
            elseif waypointEventNullCondition
                % AMEND THE META OBJECT
                SIM.OBJECTS(entityA).objectStatus(entityB,eventType.waypoint) = 0;
                % GENERATE THE WAYPOINT NULLIFICATION EVENT
                %[EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.null_waypoint);
                %metaEVENTS = vertcat(metaEVENTS,EVENT);                                        % Add new META EVENTS to global structure
            end 
            continue    % Skip collision check
        end

        % ////////// CONFIRM THE OBJECT CAN BE COLLIDED WITH //////////////      
        if SIM.OBJECTS(entityB).hitBox == OMAS_hitBoxType.none
            % Check hit-box definition:
            % -If no hitbox is assigned, no need to evaluate warnings or
            % collisions.
            continue              
        end
        
        % /////////////////// ASSESS WARNING CONDITIONS ///////////////////
        % WARNING EVENT LOGIC
        ABwarningConstraint       =  warningLogicals(entityA,entityB);
        novelWarningConstraint    =  SIM.OBJECTS(entityA).objectStatus(entityB,eventType.warning) == 0;
        warningEventCondition     =  ABwarningConstraint &&  novelWarningConstraint;         % Is within a range and a warning event has not been issued.
        warningEventNullCondition = ~ABwarningConstraint && ~novelWarningConstraint;         % Is outside a range and a warning condition is still toggled.       
        % EVALUATE OCCURANCE
        if warningEventCondition 
            % UPDATE WARNING STATUS
            SIM.OBJECTS(entityA).objectStatus(entityB,eventType.warning) = 1;
            % GENERATE THE WARNING EVENT
            [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.warning);
            metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
        elseif warningEventNullCondition
            % UPDATE WARNING-NULL STATUS
            SIM.OBJECTS(entityA).objectStatus(entityB,eventType.warning) = 0;
            % GENERATE THE WARNING NULLIFICATION EVENT
            [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.null_warning);
            metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
        end
        
        % ////////////////// ASSESS COLLISION CONDITIONS //////////////////
        % COLLISION EVENT LOGIC
        ABcollisionConstraint       =  collisionLogicals(entityA,entityB);
        novelCollisionConstraint    =  SIM.OBJECTS(entityA).objectStatus(entityB,eventType.collision) == 0;   
        collisionEventCondition     =  ABcollisionConstraint &&  novelCollisionConstraint;   % Is within a range and a collision event has not been issued.
        collisionEventNullCondition = ~ABcollisionConstraint && ~novelCollisionConstraint;   % Is outside a range and a collision condition is still toggled.     
        % EVALUATE OCCURANCE
        if collisionEventCondition
            % UPDATE COLLISIONL STATUS
            SIM.OBJECTS(entityA).objectStatus(entityB,eventType.collision) = 1;
            % GENERATE THE COLLISION EVENT
            [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.collision);
            metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
        elseif collisionEventNullCondition
        	% AMEND THE META OBJECT
            SIM.OBJECTS(entityA).objectStatus(entityB,eventType.collision) = 0;
            % GENERATE THE COLLISION NULLIFICATION EVENT
            [EVENT] = OMAS_eventHandler(SIM.TIME.currentTime,SIM.OBJECTS(entityA),SIM.OBJECTS(entityB),eventType.null_collision);
            metaEVENTS = vertcat(metaEVENTS,EVENT);                                          % Add new META EVENTS to global structure
        end     
    end
end
       
end
% UPDATE THE OBJECT PROPERTIES
function [referenceObject,objectEVENTS] = UpdateObjects(SIM,objectIndex,referenceObject)
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
SIMfirstObject = SIM.OBJECTS(SIM.globalIDvector == referenceObject.objectID);

% If the observing object is "2D"
if referenceObject.Is3D()
    dimensionIndices = 1:3;
else
    dimensionIndices = 1:2;
end

if SIM.verbosity >= 2
    fprintf('[UPDATE]\tCycling ID:%d\t name: %s\n',SIMfirstObject.objectID,referenceObject.name);

end

% Default timing parameters
ENV = SIM.TIME;
ENV.outputPath = SIM.outputPath;

if SIM.verbosity >= 2
    fprintf('[UPDATE]\tCycling ID:%d\t name: %s\n',SIMfirstObject.objectID,referenceObject.name);
end

% SWITCH BEHAVIOUR BASED ON OBJECT TYPE
switch SIMfirstObject.type
    case OMAS_objectType.agent
        % AGENT - SIMULATION/ENVIROMENTAL FEEDBACK REQUIRED %%%%%%%%%%%%%%
        % Container for agent detection packet
        observationPacket = [];             
        
        % //// PARSE THE ENVIRONMENTAL OBJECTS OBSERVABLE TO THE AGENT ////
        for ID2 = 1:SIM.totalObjects
            % GET THE SIM OBJECT OF THE EVALUATION OBJECT
            SIMsecondObject = SIM.OBJECTS(ID2);
            
            % Skip condition - Get the current status of this object
            isDetected = 1 == SIMfirstObject.objectStatus(ID2,eventType.detection);
            isSameObject = SIMfirstObject.objectID == SIMsecondObject.objectID; 
            if ~isDetected || isSameObject
                continue
            end
            
            % GET THE EQUIVALENT OBJECT TO EVALUATE AGAINST THE REFERENCE
            secondObject = objectIndex{SIM.globalIDvector == SIMsecondObject.objectID};
                            
            % ///////////////// ELSE AGENT IS DETECTED ////////////////////
            % ////// (GENERATE A LOCALISED PROXIMITY DESCIPTION) //////////
            % The second agent is within global seperation conditions. A
            % communication packet is therefore assembled containing the
            % data on the second object. rotated into the agent-local
            % coordinate frame.

            % GET THE RELATIVE POSITIONS OF THE AGENT AND OBJECT IN THE 
            % GLOBAL COORDINATE SYSTEM (x_rel = x_B - x_A)
            relativePosition = SIMfirstObject.relativePositions(ID2,:)';
            relativeVelocity = SIMsecondObject.globalState(4:6) - SIMfirstObject.globalState(4:6);
            % ROTATE THE GLOBAL STATE OF THE OBJECT INTO THE AGENT FRAME
            observedPosition = SIMfirstObject.R*relativePosition;            % Rotate the from the global into the body frame of the simReference
            observedVelocity = SIMfirstObject.R*relativeVelocity;                   
            % SPHERICAL REPRESENTATION
            observedRange     = norm(observedPosition);
            observedElevation = asin(observedPosition(3)/observedRange);
            observedHeading   = atan2(observedPosition(2),observedPosition(1));
            % PRESENT SIZE PROPERTIES
            observedRadius    = SIMfirstObject.radius;     
            observedAngularWidth = 2*asin(observedRadius/(observedRange + observedRadius));               
                        
            % OBJECT GEOMETRY PREPARATION
            % In order for the agent to observe the object correctly. The
            % geometry must be translated and rotated into the relative
            % frame of the detecting agent.
            % ASSUMPTION:
            % - The geometry is normalised to the body axes of the object.
            % TO DO:
            % - If the object is on the edge of vision.. only part of the
            %   object will be visable to the agent. a process must be in
            %   place to create a subset of the geometry.

            % THE GEOMETRIC PARAMETERS
            if size(secondObject.GEOMETRY.vertices,1) > 0
                % Rotate evaluation object geometry into the global space,
                % then rotate it back into the local space of the reference.                
                % THE RELATIVE ROTATIONS OF THE SECOND BODY
                relativeR = SIMfirstObject.R'*SIMsecondObject.R; 
                % THE CONTAINER FOR THE RELATIVE GEOMETRY
                observedGeometry = struct('vertices',secondObject.GEOMETRY.vertices*relativeR + observedPosition',...
                                          'normals', secondObject.GEOMETRY.normals*relativeR,...
                                          'faces',   secondObject.GEOMETRY.faces,...
                                          'centroid',secondObject.GEOMETRY.centroid + observedPosition');
                % [TO-DO] Reduce the structure sent to match the exposed geometry                      
                                      
%                 % PROCESS THE VISIBLE GEOMETRY (The sub-set of the geometry that is visible to the agent.)
%                 observedGeometry = OMAS_restrictedGeometry(zeros(3,1),SIMfirstObject.detectionRadius,relativeGeometry);                   
            else
                % THE OBSERVED GEOMETRY IS EMPTY
                observedGeometry = secondObject.GEOMETRY;                  % Pass the empty structure
            end
            
            % OBJECT PRIORITY (IF WAYPOINT)
            observedPriority = NaN;
            if SIMsecondObject.type == OMAS_objectType.waypoint
                association = objectIndex{ID2}.GetAgentAssociation(referenceObject);
                observedPriority = association.priority;
            end
            
            % A DEBUG STUCTURE 
            DEBUG = struct('globalPosition',SIMsecondObject.globalState(dimensionIndices),...     % Added for simplicity
                           'globalVelocity',SIMsecondObject.globalState(dimensionIndices + 3),... % Added for simplicity
                           'priority',observedPriority);                                          % The objects global priority for that agent
            
            % ///////////////// PREPARE DETECTION PACKET //////////////////
            detectionObject = struct('objectID',SIMsecondObject.objectID,...            % The object ID (observed)
                                     'name',SIMsecondObject.name,...                    % The object name tag (observed)
                                     'type',SIMsecondObject.type,...                    % The objects sim-type enum
                                     'radius',observedRadius,...                        % The objects true size
                                     'position',observedPosition(dimensionIndices,1),...% The apparent position in the relative frame
                                     'velocity',observedVelocity(dimensionIndices,1),...% The apparent velocity in the relative frame
                                     'range',observedRange,...                          % The apparent range
                                     'elevation',observedElevation,...                  % The apparent inclination angle
                                     'heading',observedHeading,...                      % The apparent Azimuth angle
                                     'width',observedAngularWidth,...                   % The apparent angular width at that range    
                                     'geometry',observedGeometry,...                    % The observable geometrical components
                                     'colour',SIMsecondObject.colour,...                % Finally the simulation's colourID 
                                     'DEBUG',DEBUG);                                    % Pass additional parameters (not otherwise known)
            observationPacket = vertcat(observationPacket,detectionObject);             % Append object to packet to agent
        end

        % /////////////// COMPUTE THE LOCAL AGENT CYCLE ///////////////////
        % Given the detection object defined for this agent, compute the
        % agents cycle with this information.
        try
            % SEND OBJECT OBSERVERATION PACKET TO AGENT
            referenceObject = referenceObject.main(ENV,observationPacket);                      % Hand the current sim TIME object and complete observation packet
        catch objectProcessError             
            warning('[%s]: Cycle failed, there was a problem with the object (%s) file, please check: %s',...
                     SIM.phase,referenceObject.name,objectProcessError.stack(1).name);
            rethrow(objectProcessError);    
        end
    case OMAS_objectType.obstacle
        % PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        referenceObject = referenceObject.main(ENV);
    case OMAS_objectType.waypoint
        % WAYPOINTS ARE CONSIDERED PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%
        referenceObject = referenceObject.main(ENV);
    otherwise
        % CONSIDERED PASSIVE - NO FEEDBACK REQUIRED %%%%%%%%%%%%%%%%%%%%%%%
        % DO NOT UPDATE
end
end

% /////////////////// SIMULATION OUTPUT OPERATIONS ////////////////////////
% GET INITIAL OUTPUT DATA STRUCTURE
function [DATA] = GetOutputStructure(SIM)
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