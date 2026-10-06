%% 2D VECTOR TOOLS (vector_tools.m) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shared planar vector operations used by avoidance algorithms.

classdef vector_tools
    methods (Static)
        function [absLength] = abs(v)
            absLength = sqrt(dot(v,v));
        end
        function [absSq] = absSq(v)
            absSq = dot(v,v);
        end
        function [det_uv] = det(u,v)
            det_uv = u(1)*v(2) - u(2)*v(1);
        end
        function [dot_uv] = dot(u,v)
            dot_uv = u(1)*v(1) + u(2)*v(2);
        end
        function [unitVector] = normalise(v)
            unitVector = v/sqrt(dot(v,v));
        end
    end
end
