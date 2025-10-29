classdef stochastic_generation < dynamicprops
    %In this class, we define methods for SGVP, SGVF, SG of Vortex cores
    %and addition. The properties of the object is included as well.
    
    properties
        u_tau, z_0, ks, delta, lambda, nu, ro_uw, N_prof, z, Delta_x,...
            Gen_u_prof, Gen_w_prof, log_data_SGVP,...
            Gen_u_LRVF, Gen_w_LRVF, Gen_x_LRVF, hist_corr,...
            Kr, Gen_u_HRVF, Gen_w_HRVF, Gen_x_HRVF,...
            Lambda_ci, Lambda_cirms,...
            ro_r_u_omega_near_wall, ro_r_u_omega_far_wall,...
            prograde_VorX_data, retrograde_VorX_data,...
            PowSpecDen_k, wavenumb, epsilon_spec, eta_spec,...
            D11, r11, epsilon_str, eta_str
    end
    
    methods
        function obj = stochastic_generation(u_tau, z_0, ks,delta, lambda, nu, ...
                z_min, z_max, res_z)
            %STOCHASTIC_GENERATION Construct an instance of this class
            obj.u_tau = u_tau;
            obj.z_0 = z_0;
            obj.ks = ks;
            obj.delta = delta;
            obj.lambda = lambda;
            obj.nu = nu;
            obj.z = z_min:res_z:z_max;
        end
        
        function SGVP(obj, ro_uw, N)
            %In this function, we generate N profiles stochastically.
            %   Detailed explanation goes here
            obj.ro_uw = ro_uw;
            obj.N_prof = N;
            [obj.Gen_u_prof, obj.Gen_w_prof, obj.log_data_SGVP] =...
                SG_VelProf(obj);

            
        end
        function SGVF(obj, Delta_x)
            %In this function, we try to reorganize the profiles generated
            %in SGVP, to have correlated velocity field.
            %   Detailed explanation goes here
            obj.Delta_x = Delta_x;
            [obj.Gen_u_LRVF, obj.Gen_w_LRVF, obj.hist_corr] =...
                SG_VelField(obj);

            
        end
        function SGVorX(obj, u_wOu_tau_near, r_wOlambda_T_near, rho_uw_VorX_near,...
            u_wOu_tau_far, r_wOlambda_T_far, rho_uw_VorX_far)
            % In this function, we try to add vortices to the SGVF.
            %   Detailed explanation goes here
            SG_VorX(obj, u_wOu_tau_near, r_wOlambda_T_near, rho_uw_VorX_near,...
                u_wOu_tau_far, r_wOlambda_T_far, rho_uw_VorX_far);
            
        end

        function Structure_function(obj)
            str_func(obj)
        end

        function Spectral_analysis(obj)
            spec_analysis(obj)    
        end
        function Gauss_kernel(obj,kern_xsize, kern_ysize, Sigma)
            % Applying Multivariate Gaussian Kernel
            MVG_kernel(obj, kern_xsize, kern_ysize, Sigma)
            
        end

    end
end

