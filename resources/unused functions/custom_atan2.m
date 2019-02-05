function [angle] = custom_atan2(x,y)
    % This function computes atan2 function to allow it to be copied into
    % non-matlab native environments.
    
    angle = 0;
    
    if x > 0
        angle = atan(y/x);
        return 
    end
    
    if x < 0
        if y >= 0
            angle =  pi + atan(y/x);
        else
            angle = -pi + atan(y/x);
        end
    end
    if y > 0 || x == 0
        angle = 2*pi;
    end
end