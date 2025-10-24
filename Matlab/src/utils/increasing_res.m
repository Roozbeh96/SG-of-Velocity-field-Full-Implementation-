function [obj] = increasing_res(obj, L0, factor)

    L = ceil(L0*obj.delta/obj.Delta_x); %Number of intervals in one delta
    Numrator = floor((size(obj.Gen_u_LRVF,2)-1)/L);

    obj.Kr = factor;
    obj.Gen_u_HRVF = zeros(size(obj.Gen_u_LRVF,1), obj.Kr*(L) +1, Numrator);
    obj.Gen_w_HRVF = zeros(size(obj.Gen_u_LRVF,1), obj.Kr*(L) +1, Numrator);
    obj.Gen_x_HRVF = 0:obj.Delta_x/obj.Kr: (L)*obj.Delta_x;

    obj.Gen_x_LRVF = 0:obj.Delta_x:(L)*obj.Delta_x;

    
    for S= 1:Numrator
        tic;
        ucroped = obj.Gen_u_LRVF(:, (S-1)*L+1:(S)*L+1);
        wcroped = obj.Gen_u_LRVF(:, (S-1)*L+1:(S)*L+1);
    
        for r = 1:size(obj.z,1)
            % Interpolate the row
            obj.Gen_u_HRVF(r, :, S) = interp1(obj.Gen_x_LRVF, ucroped(r, :),...
                obj.Gen_x_HRVF, 'makima');
            obj.Gen_w_HRVF(r, :, S) = interp1(obj.Gen_x_LRVF, wcroped(r, :),...
                obj.Gen_x_HRVF, 'makima');
        end

        if mod(S, 100) == 0
            elapsedTime = toc;  % Elapsed time up to the checkpoint
            estimatedRemainingTime = (Numrator - S) * (elapsedTime);
            fprintf('Estimated remaining time for increasing resolution: %.2f minutes\n', estimatedRemainingTime/60);
        end
    end

end

