function Lambda_ci(obj)

    obj.Lambda_ci = zeros(size(obj.Gen_u_HRVF));
    DelxGen = obj.Gen_x_HRVF(2)-obj.Gen_x_HRVF(1);
    DelzGen = obj.z(2)-obj.z(1);

    for S= 1:size(obj.Lambda_ci,3)
    
    [dudx, dudz] = gradient(obj.Gen_u_HRVF(:,:, S),DelxGen,DelzGen);
    [dwdx, dwdz] = gradient(obj.Gen_w_HRVF(:,:, S),DelxGen,DelzGen);

    omega = dwdx - dudz;
    Mat = zeros(size(dudx,1),size(dudx,2),2,2);
    Mat(:,:,1,1) = dudx;
    Mat(:,:,1,2) = dudz;
    Mat(:,:,2,1) = dwdx;
    Mat(:,:,2,2) = dwdz;
    for r = 1:size(obj.Lambda_ci,1)
    
        for c = 1:size(obj.Lambda_ci,2)
            temp = reshape(Mat(r,c,:,:),[2,2]);
            obj.Lambda_ci(r,c,S) = unique(abs(imag(eig(temp))))...
                *sign(omega(r,c));
    
        end
    end
    
    
    end
    % Since mean(Lambda_ci) ~= 0, the Lambda_cirms can be calculated as
    % follows:
    obj.Lambda_cirms = (mean(mean(obj.Lambda_ci.^2,3),2)).^0.5;
end

