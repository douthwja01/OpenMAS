% GET THE UNIT VECTOR
function [v_unit] = unit(v)
    v_unit = v/norm(v);
end