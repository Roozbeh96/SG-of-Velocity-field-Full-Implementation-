function [pr, rt] = Lambdaci_VorX(obj,i,...
    S,xcs,zcs,rot,r_wOlambda_T,...
    u_wOu_tau, rho_uw_VorX, pr, rt, N)


           
    r_omega = r_wOlambda_T*obj.lambda;
    u_omega = u_wOu_tau*obj.u_tau;
    Gama=2*pi*u_omega*r_omega/(1-exp(-1));
    Trans = [1, rho_uw_VorX; 0, sqrt(1 - rho_uw_VorX^2)];

    
    [indz] = find(obj.HRVFz>=zcs(i)-N*r_omega &...
        obj.HRVFz<=zcs(i)+N*r_omega);
    [indx] = find(obj.HRVFx>=xcs(i)-N*r_omega &...
        obj.HRVFx<=xcs(i)+N*r_omega);
    [X,Z] = meshgrid(obj.HRVFx(indx),obj.HRVFz(indz));
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
    
    mask1 = (Theta >= 0 & Theta < pi/2 & r~=0) ;
    
    mask2 =(Theta >= pi/2 & Theta <= pi+0.1 & r~=0);
    
    mask3 =((Theta >= -pi & Theta <= -pi/2 & r~=0));
    
    mask4 = (Theta >= -pi/2 & Theta <= 0 & r~=0) ;
    
    
    if rot(i)<0

        obj.prograde_VorX_data{S}{pr} = ...
        dictionary(["xc_VorX[m]","zc_VorX[m]","u_omega[m/s]","r_omega[m]",...
        "rotation[{-1,1}]"],[xcs(i), zcs(i),u_omega,r_omega,rot(i)]); 
        pr = pr +1;

        u_VorX = uazi.*sin(Theta);
        w_VorX = -uazi.*cos(Theta);

        A1 = zeros(size(w_VorX(mask1),1),2);
        A2 = zeros(size(w_VorX(mask2),1),2);
        A3 = zeros(size(w_VorX(mask3),1),2);
        A4 = zeros(size(w_VorX(mask4),1),2);
        A1(:,1) = u_VorX(mask1);
        A1(:,2) = w_VorX(mask1);
        A2(:,1) = u_VorX(mask2);
        A2(:,2) = w_VorX(mask2);
        A3(:,1) = u_VorX(mask3);
        A3(:,2) = w_VorX(mask3);
        A4(:,1) = u_VorX(mask4);
        A4(:,2) = w_VorX(mask4);
        B1 = A1* Trans;
        B2 = A2* Trans;
        B3 = A3* Trans;
        B4 = A4* Trans;
        u_TvorX = zeros(size(u_VorX));
        w_TvorX = zeros(size(w_VorX));
        u_TvorX(mask1) = B1(:,1);
        w_TvorX(mask1) = B1(:,2);     
        u_TvorX(mask2) = B2(:,1);
        w_TvorX(mask2) = B2(:,2);  
        u_TvorX(mask3) = B3(:,1);
        w_TvorX(mask3) = B3(:,2);
        u_TvorX(mask4) = B4(:,1);
        w_TvorX(mask4) = B4(:,2);
        
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




        A1 = zeros(size(w_VorX(mask1),1),2);
        A2 = zeros(size(w_VorX(mask2),1),2);
        A3 = zeros(size(w_VorX(mask3),1),2);
        A4 = zeros(size(w_VorX(mask4),1),2);
        A1(:,1) = u_VorX(mask1);
        A1(:,2) = w_VorX(mask1);
        A2(:,1) = u_VorX(mask2);
        A2(:,2) = w_VorX(mask2);
        A3(:,1) = u_VorX(mask3);
        A3(:,2) = w_VorX(mask3);
        A4(:,1) = u_VorX(mask4);
        A4(:,2) = w_VorX(mask4);
        B1 = A1* Trans;
        B2 = A2* Trans;
        B3 = A3* Trans;
        B4 = A4* Trans;
        u_TvorX = zeros(size(u_VorX));
        w_TvorX = zeros(size(w_VorX));
        u_TvorX(mask1) = B1(:,1);
        w_TvorX(mask1) = B1(:,2);
        u_TvorX(mask2) = B2(:,1);
        w_TvorX(mask2) = B2(:,2);
        u_TvorX(mask3) = B3(:,1);
        w_TvorX(mask3) = B3(:,2);
        u_TvorX(mask4) = B4(:,1);
        w_TvorX(mask4) = B4(:,2);
    



        A = [reshape(u_VorX,[],1) reshape(w_VorX,[],1)];
        B = A*TVorX;
        u_TVorX = reshape(B(:,1), size(u_VorX));
        w_TVorX = reshape(B(:,2), size(w_VorX));
    
    end
    if isempty(u_TvorX)
        return
    end
    % Adding to the velocity field
    obj.Gen_u_HRVF(indz,indx,S) = obj.Gen_u_HRVF(indz,indx,S)+ u_TvorX;
    
    obj.Gen_w_HRVF(indz,indx,S) = obj.Gen_w_HRVF(indz,indx,S)+ w_TvorX;

end

