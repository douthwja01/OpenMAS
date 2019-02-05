% THE SEPARATION RULE
function [v_sep] = separationRule(neighbours)
% We assume the position of the neighbours are measured
% relatively. We construct a vector that is the sum of
%             v_sep = zeros(3,1);
v_sep = [1;0;0];
for n = 1:numel(neighbours)
    pn = neighbours(n).position;
    % Get the (-ve) separation vector
    if sum(abs(pn)) == 0
        pn = [1;0;0]*1E-5;
    end
    vp = -pn;
    % Normalise
    vn = OMAS_geometry.unit(vp);
    % Scale by the proximity
    vn = vn*norm(vn);
    % Combine for net influence vector
    v_sep = v_sep + vn;
end

end
% THE ALIGNMENT RULE
function [v_ali] = alignmentRule(neighbours)
% This rule urges the robots to move in the same direction as
% its neighbours. (i.e heading matching)
%             v_ali = zeros(3,1);
v_ali = [1;0;0];
for n = 1:numel(neighbours)
    % THE RELATIVE VELOCITY VECTORS
    v_ali = v_ali + neighbours(n).velocity;
end
v_ali = v_ali/numel(v_ali);
end
% THE COHESION RULE
function [v_coh] = cohesionRule(neighbours)
% This rule will try to move the agent towards the center of
% mass of the other robots to form a group/swarm.
%             v_coh = zeros(3,1);
v_coh = [1;0;0];
for n = 1:numel(neighbours)
    % SUM THE VECTOR POSITIONS
    v_coh = v_coh + neighbours(n).position;
end
% THE MEAN POSITION
v_coh = v_coh/numel(neighbours);
end
% THE MIGRATION RULE
function [v_mig] = migrationRule(targetWaypoint)
% This rule aims to move the swarm/individual towards a
% designated goal location.
if ~isstruct(targetWaypoint)
    v_mig = [1;0;0];
    return
end
% THE RELATIVE POSITION VECTOR
v_mig = targetWaypoint.position;
end