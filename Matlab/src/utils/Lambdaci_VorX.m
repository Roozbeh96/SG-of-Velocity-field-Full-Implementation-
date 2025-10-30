function [pr, rt] = Lambdaci_VorX(obj,i,...
    S,xcs,zcs,rot,r_wOlambda_T,...
    u_wOu_tau, rho_uw_VorX, pr, rt, N)


           
    r_omega = r_wOlambda_T*obj.lambda;
    u_omega = u_wOu_tau*obj.u_tau;
    Gama=2*pi*u_omega*r_omega/(1-exp(-1));
    TVorX = [1, rho_uw_VorX; 0, sqrt(1 - rho_uw_VorX^2)];

    
    [indz] = find(obj.z>=zcs(i)-N*r_omega &...
        obj.z<=zcs(i)+N*r_omega);
    [indx] = find(obj.Gen_x_HRVF>=xcs(i)-N*r_omega &...
        obj.Gen_x_HRVF<=xcs(i)+N*r_omega);
    [X,Z] = meshgrid(obj.Gen_x_HRVF(indx),obj.z(indz));
    r = sqrt((X-xcs(i)).^2+(Z-zcs(i)).^2);
    % Making the vortex circular shape.
    r(r>N*r_omega)=0;
    
    if size(r,1)==1
        return
    end
    
    
    uazi=Gama./(2*pi*r).*(1-exp(-r.^2./r_omega^2));
    nan_logical_array = isnan(uazi);
    [nan_indices] = find(nan_logical_array);
    uazi(nan_indices) = 0;
    
    Theta = atan2((Z-zcs(i)),(X-xcs(i)));
    nan_logical_array = isnan(Theta);
    [nan_indices] = find(nan_logical_array);
    Theta(nan_indices) = 0;
    
    
    if rot(i)<0

        obj.prograde_VorX_data{S}{pr} = ...
        dictionary(["xc_VorX[m]","zc_VorX[m]","u_omega[m/s]","r_omega[m]",...
        "rotation[{-1,1}]"],[xcs(i), zcs(i),u_omega,r_omega,rot(i)]); 
        pr = pr +1;

        u_VorX = uazi.*sin(Theta);
        w_VorX = -uazi.*cos(Theta);

        if isempty(u_VorX)
            return
        end
        
        
        A = [reshape(u_VorX,[],1) reshape(w_VorX,[],1)];
        B = A*TVorX;
        u_TVorX = reshape(B(:,1), size(u_VorX));
        w_TVorX = reshape(B(:,2), size(w_VorX));
        
    else
      
        obj.retrograde_VorX_data{S}{rt} = ...
        dictionary(["xc_VorX[m]","zc_VorX[m]","u_omega[m/s]","r_omega[m]",...
        "rotation[{-1,1}]"],[xcs(i), zcs(i),u_omega,r_omega,rot(i)]); 
        rt = rt+1;
        u_VorX = -uazi.*sin(Theta);
        w_VorX = uazi.*cos(Theta);

        if isempty(u_VorX)
            return
        end

        A = [reshape(u_VorX,[],1) reshape(w_VorX,[],1)];
        B = A*TVorX;
        u_TVorX = reshape(B(:,1), size(u_VorX));
        w_TVorX = reshape(B(:,2), size(w_VorX));
    
    end
    
    % Adding to the velocity field
    obj.Gen_u_HRVF(indz,indx,S) = obj.Gen_u_HRVF(indz,indx,S)+ u_TVorX;
    
    obj.Gen_w_HRVF(indz,indx,S) = obj.Gen_w_HRVF(indz,indx,S)+ w_TVorX;

end

