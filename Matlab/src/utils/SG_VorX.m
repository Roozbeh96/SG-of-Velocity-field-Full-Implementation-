function [obj] = SG_VorX(obj)
    % At the first, step we need to divide the VF into 1 * delta (L = 1)
    % intervals, and increase the resolution. For increasing the resolution,
    % the factor that should be increased is another input parameter. Here 
    % Kr is the factor parameter which is consistent with the parameter in 
    % the paper.  All these steps are done increasing_res function.
    L = 1; 
    Kr = 32;
    [obj] = increasing_res(obj, L, Kr);

end

