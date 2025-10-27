function Lambda_ci(obj)

    obj.Lambda_ci = zeros(size(obj.Gen_u_HRVF));
    DelxGen = obj.Gen_x_HRVF(2)-obj.Gen_x_HRVF(1);
    DelzGen = obj.z(2)-obj.z(1);

    for S= 1:size(obj.Lambda_ci,3)
    tic;
    [dudx, dudz] = gradient(obj.Gen_u_HRVF(:,:, S),DelxGen,DelzGen);
    [dwdx, dwdz] = gradient(obj.Gen_w_HRVF(:,:, S),DelxGen,DelzGen);

    omega = dwdx - dudz;
        for r = 1:size(obj.Lambda_ci,1)
        
            for c = 1:size(obj.Lambda_ci,2)
    
                temp = [dudx(r,c) dudz(r,c); dwdx(r,c) dwdz(r,c)];
                obj.Lambda_ci(r,c,S) = unique(abs(imag(eig(temp))))...
                    *sign(omega(r,c));
        
            end
        end
        if mod(S, 100) == 0
                elapsedTime = toc;  % Elapsed time up to the checkpoint
                estimatedRemainingTime = (size(obj.Lambda_ci,3) - S) * (elapsedTime);
                fprintf('Estimated remaining time for computing Lambda_ci field: %.2f minutes\n', estimatedRemainingTime/60);
        end
      
    end
    
    obj.Lambda_cirms = mean(std(obj.Lambda_ci,0,3),2);
end

