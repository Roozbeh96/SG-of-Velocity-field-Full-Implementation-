function SG_VorX(obj, u_wOu_tau_near, r_wOlambda_T_near,...
    u_wOu_tau_far, r_wOlambda_T_far)
        %{
            The Pearson correlation coefficient between radius and the
            maximum azimuthal velocity of the generated vortices at
            different wall normal elevations.
            Near : z/k_s < 2
            Far : z/k_s >= 2
        %}
        rho_near = 0.4; 
        rho_Far = 0.45; 


        ii = 1;
        iii = 1;
        z_prog = zeros(1,1);
        z_retro = zeros(1,1);
        z_prog_near = zeros(1,1);
        pr = 1;
        rt = 1;
        pr_near = 1;
        rho_Near_wall = load('samples_Near.mat').samples_Near;
        rho_Far_wall = load('samples_Far.mat').samples_Far;
        xi = 1;
        xii = 1;
        Dist_close_wall_VorX = zeros(1);
        ind_close_wall_VorX = 1;
        Lambda_ci_temp = obj.Lambda_ci;
        for r = 1:size(obj.Lambda_ci, 1)
                Lambda_ci_temp(r, :, :) = obj.Lambda_ci(r, :, :) .* (abs(obj.Lambda_ci(r, :, :)) >= 1.0*obj.Lambda_cirms(r));
        end
        if obj.name == "GenASL"
            x_ = obj.Delx; %#ok<NASGU>
            obj.Delx = 0.15;
        end
        for S=1:size(obj.HRVFu,3)
            xend_near_wall_prev = 0;
            temp_ind_close_wall_VorX = 1;
            while xend_near_wall_prev < max(obj.HRVFx)% Attached wall vortices

                D1 = 1;
                D2 = 1;
                N = 1;
                fac = 1.0;% should be more than 1.0
                if obj.name ~= "GenASL"
                    r_omega = r_wOlambda_T_near(ii)*obj.Delx; %choose from near wall distribution
%                         r_omega = rwOlambda_T_near(ii)*mean(obj.Lambda_T); %choose from near wall distribution
                    Trans = [1, rho_Near_wall(xi); 0, sqrt(1 - rho_Near_wall(xi)^2)];
                    Gama=2*pi*u_wOu_tau_near(ii)*obj.u_tau*r_omega/(1-exp(-1)); 
                    xi= xi+1;
                    ii = ii+1;
                else
                    r_omega = r_wOlambda_T_far(iii)*obj.Delx;
%                         r_omega = rwOlambda_T_Far(iii)*obj.Lambda_T;
                    Trans = [1, rho_Far_wall(xii); 0, sqrt(1 - rho_Far_wall(xii)^2)];
                    Gama=2*pi*u_wOu_tau_far(iii)*obj.u_tau*r_omega/(1-exp(-1)); 
                    xii = xii+1;
                    iii = iii+1;
                end
                zc_near_wall = N*r_omega + obj.HRVFz(1);
                xc_near_wall = N*r_omega + xend_near_wall_prev;
                xend_near_wall_prev = xc_near_wall + fac*N*r_omega;
                [indzmid] = find(obj.HRVFz>=zc_near_wall,1);
                [indxmid] = find(obj.HRVFx>=xc_near_wall,1);


                %                     logic = randi([0, 1]);

                if obj.HRVFuprime(indzmid,indxmid,S)<0
%                         ii = ii + 1;
%                         xi= xi+1;
%                         xii = xii+1;
                    continue
                end

                if randi([0, 1])==0
%                         ii = ii + 1;
%                         xi= xi+1;
%                         xii = xii+1;
                    continue
                end

                if temp_ind_close_wall_VorX > 1
                    Dist_close_wall_VorX(ind_close_wall_VorX-1) = xc_near_wall-prev_pos_xc_near_wall_VorX;
                end
                prev_pos_xc_near_wall_VorX = xc_near_wall;
                temp_ind_close_wall_VorX = temp_ind_close_wall_VorX +1;
                ind_close_wall_VorX = ind_close_wall_VorX +1;

                %                     if zc_near_wall<= 2*obj.ks

                %                     Trans = [1, rho_Near_wall(xi); 0, sqrt(1 - rho_Near_wall(xi)^2)];
                %                     xi= xi+1;
                %                     else

                %                         Trans = [1, rho_Far_wall(xii); 0, sqrt(1 - rho_Far_wall(xii)^2)];
                %                         xii = xii+1;
                %                     end

%                     Gama=2*pi*delUwOu_tau_near(ii)*obj.u_tau*r_omega/(1-exp(-1));  %choose from near wall distribution
%                     ii = ii + 1;
                [indz] = find(obj.HRVFz>=zc_near_wall-N*r_omega &...
                    obj.HRVFz<=zc_near_wall+N*r_omega);
                [indx] = find(obj.HRVFx>=xc_near_wall-N*r_omega &...
                    obj.HRVFx<=xc_near_wall+N*r_omega);
                [X,Z] = meshgrid(obj.HRVFx(indx),obj.HRVFz(indz));
                r = sqrt((X-xc_near_wall).^2+(Z-zc_near_wall).^2);
                %                     r(r>N*1.01*r_omega)=0;
                z_prog_near(pr_near,1) = zc_near_wall;

                pr_near = pr_near +1;

                uazi=Gama./(2*pi*r).*(1-exp(-r.^2./r_omega^2));
                nan_logical_array = isnan(uazi);
                [nan_indices] = find(nan_logical_array);
                uazi(nan_indices) = 0;

                Theta = atan2((Z-zc_near_wall),(X-xc_near_wall));
                nan_logical_array = isnan(Theta);
                [nan_indices] = find(nan_logical_array);
                Theta(nan_indices) = 0;

                mask1 = (Theta >= 0 & Theta < pi/2 & r~=0) ;

                mask2 =(Theta >= pi/2 & Theta <= pi+0.1 & r~=0);

                mask3 =((Theta >= -pi & Theta <= -pi/2 & r~=0));

                mask4 = (Theta >= -pi/2 & Theta <= 0 & r~=0) ;


                u_vor = uazi.*sin(Theta);
                u_vor(uazi~=0)=u_vor(uazi~=0)-Gama/(2*pi*r_omega)*(1-exp(-1));
                w_vor = -uazi.*cos(Theta);
                %                     w_vor(uazi>0)=w_vor(uazi>0)-Gama/(2*pi*r_omega)*(1-exp(-1));
                A1 = zeros(size(w_vor(mask1),1),2);
                A2 = zeros(size(w_vor(mask2),1),2);
                A3 = zeros(size(w_vor(mask3),1),2);
                A4 = zeros(size(w_vor(mask4),1),2);
                A1(:,1) = u_vor(mask1);
                A1(:,2) = w_vor(mask1);
                A2(:,1) = u_vor(mask2);
                A2(:,2) = w_vor(mask2);
                A3(:,1) = u_vor(mask3);
                A3(:,2) = w_vor(mask3);
                A4(:,1) = u_vor(mask4);
                A4(:,2) = w_vor(mask4);


                B1 = A1* Trans;
                B2 = A2* Trans;
                B3 = A3* Trans;
                B4 = A4* Trans;
                u_vor_trans = zeros(size(u_vor));
                w_vor_trans = zeros(size(w_vor));
                u_vor_trans(mask1) = B1(:,1);
                w_vor_trans(mask1) = B1(:,2);

                u_vor_trans(mask2) = B2(:,1);
                w_vor_trans(mask2) = B2(:,2);

                u_vor_trans(mask3) = B3(:,1);
                w_vor_trans(mask3) = B3(:,2);

                u_vor_trans(mask4) = B4(:,1);
                w_vor_trans(mask4) = B4(:,2);

                %                     figure
                %                     plot(B1(:,1)/(delUwOu_tau_near(ii)*obj.u_tau),B1(:,2)/(delUwOu_tau_near(ii)*obj.u_tau),'r.')
                %                     hold on
                %                     plot(B2(:,1)/(delUwOu_tau_near(ii)*obj.u_tau),B2(:,2)/(delUwOu_tau_near(ii)*obj.u_tau),'r.')
                %                     plot(B3(:,1)/(delUwOu_tau_near(ii)*obj.u_tau),B3(:,2)/(delUwOu_tau_near(ii)*obj.u_tau),'r.')
                %                     plot(B4(:,1)/(delUwOu_tau_near(ii)*obj.u_tau),B4(:,2)/(delUwOu_tau_near(ii)*obj.u_tau),'r.')
                %                     set(gca,'TickLabelInterpreter','latex','FontSize',13,...
                %                         'XGrid','on','YGrid','on')
                %                     xlabel('u/$u_{\omega}$','Interpreter','Latex','FontSize',14);
                %                     ylabel('w/$u_{\omega}$','Interpreter','Latex','FontSize',14);
                %                     axis equal
                %                     xlim([-1.1 1.1])
                %                     ylim([-1.1 1.1])
                %
                %                     figure
                %                     quiver((X(1:2:end,1:2:end)-xc_near_wall)./r_omega,...
                %                         (Z(1:2:end,1:2:end)-zc_near_wall)./r_omega,...
                %                         u_vor_trans(1:2:end,1:2:end),w_vor_trans(1:2:end,1:2:end)...
                %                         ,1,'color',[0.86,0.09,0.24],...
                %                         'MaxHeadSize',0.2,'LineWidth',1.5);
                %                     set(gca,'TickLabelInterpreter','latex','FontSize',13,...
                %                         'XGrid','on','YGrid','on')
                %                     xlabel('x/r$_{\omega}$','Interpreter','Latex','FontSize',14);
                %                     ylabel('z/r$_{\omega}$','Interpreter','Latex','FontSize',14);
                %                     axis equal
                %                     xlim([-1.2 1.2])
                %                     ylim([-1.2 1.2])
                if isempty(u_vor)
                    continue
                end

                obj.HRVFu(indz,indx,S) = obj.HRVFu(indz,indx,S)+ D1*u_vor_trans;
                %                     obj.HRVFuprime(indz,indx,S) = obj.HRVFuprime(indz,indx,S)+ u_vor_trans;
                obj.HRVFw(indz,indx,S) = obj.HRVFw(indz,indx,S)+ D2*w_vor_trans;


                %                     obj.HRVFu(indz,indx,S) = obj.HRVFu(indz,indx,S)+ D1*u_vor;
                % %                     obj.HRVFuprime(indz,indx,S) = obj.HRVFuprime(indz,indx,S)+ u_vor;
                %                     obj.HRVFw(indz,indx,S) = obj.HRVFw(indz,indx,S)+ D2*w_vor;
            end
            % Lambda_ci vortices
            %Location of the vortex now is at the center(it can be stochastic
            %or based onthe swirling motion)
            %Number of the vortices should be changed based on the
            %number of lambda_ci
            Matrix = Lambda_ci_temp(:,:,S);

            binaryMatrix = Matrix ~= 0;

            % Label connected components
            [labeledMatrix, numComponents] = bwlabel(binaryMatrix, 4);  % 4-connectivity

            % Initialize variables to store centroids
            centroids = zeros(1,5);
            ix = 1;
            % Calculate the mass center (centroid) for each connected component
            for i = 1:numComponents
                % Extract the indices of the current component
                [rows, cols] = find(labeledMatrix == i);
                values = Matrix(labeledMatrix == i);
                sign_region_mean = sign(mean(values));
                % Calculate the centroid
                row_centroid = sum(rows .* values) / sum(values);
                col_centroid = sum(cols .* values) / sum(values);

                % Round the centroids to the nearest integer
                row_centroid = round(row_centroid);
                col_centroid = round(col_centroid);
                if row_centroid <= 0 || row_centroid > size(obj.HRVFz,1)
                    continue
                end
                if col_centroid <= 0 || col_centroid > size(obj.HRVFx,2)
                    continue
                end
                % Store the centroid
                centroids(ix, :) = [row_centroid, col_centroid, sign_region_mean,...
                    min(cols), max(cols)];
                ix = ix + 1;
            end

            xcs = obj.HRVFx(centroids(:, 2));

            zcs = obj.HRVFz(centroids(:, 1));

            rot = centroids(:, 3);
            Numbvortices  = size(centroids,1);
            %                 indxcs = randi([1,size(obj.HRVFu,2)],Numbvortices,1);
            %                 xcs = obj.HRVFx(indxcs);
            %                 indzcs = randi([1,size(obj.HRVFu,1)],Numbvortices,1);
            %                 zcs = obj.HRVFz(indzcs);
            %                 xcs = obj.HRVFx(floor(size(obj.HRVFu,2)/2));
            %                 zcs = obj.HRVFz(floor(size(obj.HRVFu,1)/2));
            %                 kappa = 0.39;
            for i=1:Numbvortices
                if zcs(i)<= 2*obj.ks %#ok<IFBDUP>
                    N1 = 2;
                else
                    N1 = 2;
                end
                %2This should be 2 otherwise E1, E2 has no significant effect
                E1 = 1;
                E2 = 1;
                [obj,ii,iii,z_prog, z_retro, pr, rt, xi, xii] = ...
                    lambdacivorX(obj,i,ii,iii,S,xcs,zcs,rot,r_wOlambda_T_near,r_wOlambda_T_far,...
                    u_wOu_tau_near,u_wOu_tau_far, z_prog, z_retro, pr, rt, N1, E1, E2,rho_Near_wall(xi),xi,...
                    rho_Far_wall(xii), xii, rho_near, rho_Far);%Primary vortex
                xold = xcs(i);
                if zcs(i)<= 2*obj.ks %#ok<IFBDUP>
                    N2 = 2;
                else
                    N2 = 2;
                end
                %filling periphery of the vortex(Secondary vortices)
                if zcs(i)<= 2*obj.ks % close to the wall
                    while  xold+N1*r_wOlambda_T_near(ii-1)*obj.Delx+N2*r_wOlambda_T_near(ii)*obj.Delx<obj.HRVFx(centroids(i, 5))
                        [obj,ii,iii,z_prog, z_retro, pr, rt, xi, xii] = ...
                            lambdacivorX(obj,i,ii,iii,S,xcs,zcs,rot,r_wOlambda_T_near,r_wOlambda_T_far,...
                            u_wOu_tau_near,u_wOu_tau_far, z_prog, z_retro, pr, rt, N1, E1, E2,rho_Near_wall(xi),xi,...
                            rho_Far_wall(xii), xii, rho_near, rho_Far);
                        xold = xold+N1*r_wOlambda_T_near(ii-1)*obj.Delx+N2*r_wOlambda_T_near(ii)*obj.Delx;
                        N1 = N2;
                    end
                    xold = xcs(i);
                    N1 = 2;%1
                    ii = ii +1;
                    while xold-N1*r_wOlambda_T_near(ii-1)*obj.Delx-N2*r_wOlambda_T_near(ii)*obj.Delx>obj.HRVFx(centroids(i, 4))
                        [obj,ii,iii,z_prog, z_retro, pr, rt, xi, xii] = ...
                            lambdacivorX(obj,i,ii,iii,S,xcs,zcs,rot,r_wOlambda_T_near,r_wOlambda_T_far,...
                            u_wOu_tau_near,u_wOu_tau_far, z_prog, z_retro, pr, rt, N1, E1, E2,rho_Near_wall(xi),xi,...
                            rho_Far_wall(xii), xii,rho_near, rho_Far);
                        xold = xold-N1*r_wOlambda_T_near(ii-1)*obj.Delx-N2*r_wOlambda_T_near(ii)*obj.Delx;
                        N1 = N2;
                    end
                else %Far from the wall
                    while  xold+N1*r_wOlambda_T_far(iii-1)*obj.Delx+N2*r_wOlambda_T_far(iii)*obj.Delx<obj.HRVFx(centroids(i, 5))
                        [obj,ii,iii,z_prog, z_retro, pr, rt, xi, xii] = ...
                            lambdacivorX(obj,i,ii,iii,S,xcs,zcs,rot,r_wOlambda_T_near,r_wOlambda_T_far,...
                            u_wOu_tau_near,u_wOu_tau_far, z_prog, z_retro, pr, rt, N1, E1, E2,rho_Near_wall(xi),xi,...
                            rho_Far_wall(xii), xii,rho_near, rho_Far);
                        xold = xold+N1*r_wOlambda_T_far(iii-1)*obj.Delx+N2*r_wOlambda_T_far(iii)*obj.Delx;
                        N1 = N2;
                    end
                    xold = xcs(i);
                    N1 = 2;%1
                    iii = iii +1;
                    while xold-N1*r_wOlambda_T_far(iii-1)*obj.Delx-N2*r_wOlambda_T_far(iii)*obj.Delx>obj.HRVFx(centroids(i, 4))
                        [obj,ii,iii,z_prog, z_retro, pr, rt, xi, xii] = ...
                            lambdacivorX(obj,i,ii,iii,S,xcs,zcs,rot,r_wOlambda_T_near,r_wOlambda_T_far,...
                            u_wOu_tau_near,u_wOu_tau_far, z_prog, z_retro, pr, rt, N1, E1, E2,rho_Near_wall(xi),xi,...
                            rho_Far_wall(xii), xii,rho_near, rho_Far);
                        xold = xold-N1*r_wOlambda_T_far(iii-1)*obj.Delx-N2*r_wOlambda_T_far(iii)*obj.Delx;
                        N1 = N2;
                    end
                end
            end
        end

end

