function [obj] = increasing_res(obj, L0, factor)

    L = ceil(L0*obj.delta/obj.Delta_x); %Number of intervals in one delta
    Numrator = floor((size(obj.Gen_u_LRVF,2)-1)/L);

    obj.Kr = factor;
    obj.Gen_u_HRVF = zeros(size(obj.u,1), obj.Kr*(L) +1, Numrator);
    obj.Gen_w_HRVF = zeros(size(obj.u,1), obj.Kr*(L) +1, Numrator);
    obj.HRVFx = 0:obj.Delta_x/obj.Kr: (L)*obj.Delta_x;
    obj.HRVFz = obj.z;
    xorg = 0:obj.Delta_x:(L)*obj.Delta_x;
    obj.x = xorg;
    xi = obj.HRVFx;
    
    for S= 1:Numrator
    
        ucroped = obj.u(:, (S-1)*L+1:(S)*L+1);
        wcroped = obj.w(:, (S-1)*L+1:(S)*L+1);
    
        for r = 1:size(obj.HRVFz,1)
            % Interpolate the row
            obj.Gen_u_HRVF(r, :, S) = interp1(xorg, ucroped(r, :), xi, 'makima');
            obj.Gen_w_HRVF(r, :, S) = interp1(xorg, wcroped(r, :), xi, 'makima');
        end
    end

end

