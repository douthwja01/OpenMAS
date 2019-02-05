%                 for i = 1:numel(obstacleNeighbours)
%                     % GET THE OBSTACLE
% 				    obstacle1 = obstacleNeighbours(i);
% 				    obstacle2 = obstacle1.nextObstacle;
%                     relativePosition1 = obstacle1.point - agentPosition;
%                     relativePosition2 = obstacle2.point - agentPosition;
%                     
%                     % Check if velocity obstacle of obstacle is already taken care of by
%                     % previously constructed obstacle ORCA lines.
%                     %
%                     alreadyCovered = logical(false);
% 
%                     for j = 1:numel(obj.orcaLines)
%                         conditionA = obj.det(invTimeHorizonObst*relativePosition1 - obj.orcaLines(j).point, obj.orcaLines(j).direction) - invTimeHorizonObst*agentRadius >= -obj.RVO_EPSILON;
%                         conditionB = obj.det(invTimeHorizonObst*relativePosition2 - obj.orcaLines(j).point, obj.orcaLines(j).direction) - invTimeHorizonObst*agentRadius >= -obj.RVO_EPSILON;
%                         if conditionA && conditionB
%                             alreadyCovered = logical(true);
%                             break;
%                         end
%                     end
%                     % If covered, move to the next obstacle
%                     if alreadyCovered
%                         continue
%                     end
%                     
%                     % Obstacle is not yet covered, check for collisions
%                     distSq1 = obj.absSq(relativePosition1);
%                     distSq2 = obj.absSq(relativePosition2);
%                     radiusSq = agentRadius^2;
%                     
%                     obstacleVector = obstacle2.point - obstacle1.point;
%                     s = dot(-relativePosition1,obstacleVector)/obj.absSq(obstacleVector);
%                     distSqlineObj = obj.absSq(-relativePosition1 - s*obstacleVector);
%                     
%                     % CREATE ORCA lineObj
%                     lineObj = struct('point',[],'direction',[]);
%                     
%                     if s < 0 && distSq1 <= radiusSq
%                         % Collision with left vertex. Ignore if non-convex
%                         if obstacle1.isConvex
%                             lineObj.point = [0;0];
%                             direction = [-relativePosition1(2);relativePosition1(1)];  % [y;x]
%                             lineObj.direction = obj.normalise(direction);
%                             obj.orcaLines = [obj.orcaLines;lineObj];
%                         end
%                         continue
%                     elseif s > 1 && distSq2 <=  radiusSq
%                         % Collision with right vertex. Ignore if non-convex or
%                         % if it will be taken care of by neighbouring obstacle.
%                         if obstacle2.isConvex && obj.det(relativePosition2,obstacle2.unitDir) >= 0
%                             lineObj.point = [0;0];
%                             direction = [-relativePosition2(2);relativePosition2(1)];
%                             lineObj.direction = obj.normalise(direction);
%                             obj.orcaLines = [obj.orcaLines;lineObj];
%                         end
%                         continue
%                     elseif s >= 0 && s < 0 && distSqlineObj <= radiusSq
%                         % Collision with obstacle segment
%                         lineObj.point = [0;0];
%                         lineObj.direction = -obstacle1.unitDir;
%                         obj.orcaLines = [obj.orcaLines;lineObj];
%                         continue
%                     end
%                     
%                     % NO COLLISION
%                     % Comput legs. When obliquely viewed, both legs can come
%                     % from a single vertex. Legs extend cut-off lineObj when
%                     % non-convex vertex
%                     
%                     if s < 0 && distSqlineObj <= radiusSq
%                         % Obstacle viewed obliquely so that left vertex defines
%                         % velocity obstacle
%                         if ~obstacle1.isConvex
%                             % Ignore obstacle
%                             continue;
%                         end
%                         obstacle2 = obstacle1;
%                         
%                         leg1 = sqrt(distSq1 -radiusSq);
%                         leftLegDirection  = [relativePosition1(1)*leg1 - relativePosition1(2)*obj.radius;  relativePosition1(1)*obj.radius + relativePosition1(2)*leg1] / distSq1;
%                         rightLegDirection = [relativePosition1(1)*leg1 + relativePosition1(2)*obj.radius; -relativePosition1(1)*obj.radius + relativePosition1(2)*leg1] / distSq1; % Assuming [x;y]
%                     elseif s > 1 && distSqlineObj <= radiusSq
%                         % Obstacle viewed obliquely so that right vertex defines
%                         % velocity obstacle.
%                         if ~obstacle2.isConvex
%                             % Ignore obstacle
%                             continue
%                         end
%                         obstacle1 = obstacle2;
%                         
%                         leg2 = sqrt(distSq2 - radiusSq);
%                         leftLegDirection  = [relativePosition2(1)*leg2 - relativePosition2(2)*obj.radius;  relativePosition2(1)*obj.radius + relativePosition2(2)*leg2] / distSq2;
%                         rightLegDirection = [relativePosition2(1)*leg2 + relativePosition2(2)*obj.radius; -relativePosition2(1)*obj.radius + relativePosition2(2)*leg2] / distSq2; % Assuming [x;y]
%                     else
%                         % Unusual situation
%                         if obstacle1.isConvex
%                             leg1 = sqrt(distSq1 - radiusSq);
%                             leftLegDirection = [relativePosition1(1)*leg1 - relativePosition1(2)*obj.radius; relativePosition1(1)*obj.radius + relativePosition1(2)*leg1] / distSq1;
%                         else
%                             % Left vertex non-convex; left leg extends cut-off
%                             % lineObj
%                             leftLegDirection = -obstacle1.unitDir;
%                         end
%                         
%                         if obstacle2.isConvex
%                             leg2 = sqrt(distSq2 - radiusSq);
%                             rightLegDirection = [relativePosition2(1)*leg2 + relativePosition2(2)*obj.radius; -relativePosition2(1)*obj.radius + relativePosition2(2)*leg2] / distSq2;
%                         else
%                             % Right vertex non-convex; right leg extends cut-off
%                             % lineObj
%                             rightLegDirection = obstacle1.unitDir;
%                         end
%                     end
%                     
%                     % Legs can never point into neighbouring edge when convex
%                     % vertex, take cutoff-lineObj of neighbouring edge instead. If
%                     % velocity projected on "foreign" leg, no constraint is
%                     % added.
%                     
%                     leftNeighbour = obstacle1.prevObstacle; 
%                     isLeftLegForeign  = logical(false);
%                     isRightLegForeign = logical(false);
%                     
%                     if obstacle1.isConvex && obj.det(leftLegDirection,-leftNeighbour.unitDir) >= 0
%                         % Left leg points into obstacle.
%                         leftLegDirection = -leftNeighbour.unitDir;
%                         isLeftLegForeign = logical(true);
%                     end
%                     if obstacle2.isConvex && obj.det(rightLegDirection, obstacle2.unitDir) <= 0
%                         % Right leg points into obstacle.
%                         rightLegDirection = obstacle2.unitDir;
%                         isRightLegForeign = logical(true);
%                     end
%                     % COMPUTE CUT-OFF CENTERS
%                     leftCutoff  = invTimeHorizonObst * (obstacle1.point - obj.position);
%                     rightCutoff = invTimeHorizonObst * (obstacle2.point - obj.position);
%                     cutoffVec   = rightCutoff - leftCutoff;
%                     
%                     % Project current velocity on velocity obstacle.
%                     
%                     % Check if current velocity is projected on cutoff circles
%                     if obstacle1.id == obstacle2.id
%                         t = 0.5;
%                     else
%                         t = dot((obj.velocity - leftCutoff),cutoffVec) / obj.absSq(cutoffVec);
%                     end
%                     
%                     tLeft = dot((obj.velocity - leftCutoff),leftLegDirection);
%                     tRight = dot((obj.velocity - rightCutoff),rightLegDirection);
%                     
%                     if ((t < 0 && tLeft < 0) || (obstacle1.id == obstacle2.id && tLeft < 0 && tRight < 0))
%                         % Project on left cut-off circle.
%                         unitW = obj.normalise(obj.velocity - leftCutoff);
%                         
%                         lineObj.direction = [unitW(2), -unitW(1)];
%                         lineObj.point = leftCutoff + obj.radius*invTimeHorizonObst*unitW;
%                         obj.orcaLines = [obj.orcaLines;lineObj];
%                         continue
%                     elseif (t > 1 && tRight < 0)
%                         % Project on right cut-off circle.
%                         unitW = obj.normalise(obj.velocity - rightCutoff);
%                         
%                         lineObj.direction = [unitW(2), -unitW(1)];
%                         lineObj.point = rightCutoff + agentRadius*invTimeHorizonObst*unitW;
%                         obj.orcaLines = [obj.orcaLines;lineObj];
%                         continue;
%                     end
%                     
%                     % Project on left leg, right leg, or cut-off lineObj, whichever is closest
%                     % to velocity.
%                     % CALCULATE distSqCutoff
%                     
%                     %distSqCutoff = ((t < 0.0f || t > 1.0f || obstacle1 == obstacle2) ? std::numeric_limits<float>::infinity() : absSq(velocity_ - (leftCutoff + t * cutoffVec)));
%                     conditionA = (t < 0 || t > 1 || obstacle1.id == obstacle2.id);
%                     if conditionA == 1
%                         distSqCutOff = inf;
%                     else
%                         distSqCutOff = obj.absSq(obj.velocity - (leftCutoff + t * cutoffVec));
%                     end
%                     
%                     % CALCULATE distSqLeft
%                     if (tLeft < 0) 
%                         distSqLeft = inf;
%                     else
%                         distSqLeft = obj.absSq(obj.velocity - (leftCutoff + tLeft*leftLegDirection));
%                     end
%                     % CALCULATE distSqRight
%                     if (tRight < 0)
%                         distSqRight = inf;
%                     else
%                         distSqRight = obj.absSq(obj.velocity - (rightCutoff + tRight*rightLegDirection));
%                     end 
%                     
%                     if distSqCutOff <= distSqLeft && distSqCutOff <= distSqRight
%                         % Project on cut-off lineObj.
%                         lineObj.direction = -obstacle1.unitDir;
%                         lineObj.point = leftCutoff + agentRadius*invTimeHorizonObst*[-lineObj.direction(2), lineObj.direction(1)];
%                         obj.orcaLines = [obj.orcaLines;lineObj];
%                         continue;
%                     elseif distSqLeft <= distSqRight
%                         % Project on left leg.
%                         if isLeftLegForeign
%                             continue
%                         end
%                         lineObj.direction = leftLegDirection;
%                         lineObj.point = leftCutoff + agentRadius*invTimeHorizonObst*[-lineObj.direction(2), lineObj.direction(1)];
%                         obj.orcaLines = [obj.orcaLines;lineObj];
%                         continue;
%                     else
%                         % Project on right leg. */
%                         if isRightLegForeign
%                             continue;
%                         end
%                         
%                         lineObj.direction = -rightLegDirection;
%                         lineObj.point = rightCutoff + agentRadius*invTimeHorizonObst*[-lineObj.direction(2), lineObj.direction(1)];
%                         obj.orcaLines = [obj.orcaLines;lineObj];
%                         continue;
%                     end
%                 end