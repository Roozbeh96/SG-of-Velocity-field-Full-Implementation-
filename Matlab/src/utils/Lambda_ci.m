function Lambda_ci(obj)

    obj.Lambda_ci = zeros(size(obj.Gen_u_HRVF));
    DelxGen = obj.Gen_x_HRVF(2)-obj.Gen_x_HRVF(1);
    DelzGen = obj.z(2)-obj.z(1);
    progressbar('Computing Lambda_{ci} Field')
    for S= 1:size(obj.Lambda_ci,3)
    
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
        progressbar((S)/size(obj.Lambda_ci,3))
      
    end
    
    obj.Lambda_cirms = mean(std(obj.Lambda_ci,0,3),2);
end

