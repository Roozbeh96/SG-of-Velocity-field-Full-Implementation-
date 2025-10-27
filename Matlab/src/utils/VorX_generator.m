function [u_wOu_tau_near, r_wOlambda_T_near,...
    u_wOu_tau_far, r_wOlambda_T_far]= VorX_generator(obj)

    % Generating corrolated random numbers for radius and max-azimuthal
    % velocity. First we generate for near wall then generate for far from
    % the wall.

    obj.ro_r_u_omega_near_wall = 0.4;
    n = 1e7; % Number of samples
    repro_numb = 97;
    [rand_u_omega, rand_r_omega] = genGaussCop(obj.ro_r_u_omega_near_wall,...
        n, repro_numb);
    %{
        Parameters can be found in Table 4 of the paper:
        Stochastic generation of velocity fields to reproduce the energy
        spectrum in wall turbulence
    %}
    % Generating u_omega
    a = 0.65;
    b = 4.5;
    m = 0.55;
    u_wOu_tau_near = -m*log(-rand_u_omega*(exp(-a/m)-exp(-b/m))+exp(-a/m));

    % Generating r_omega
    mu = -1.55;
    sigma = 0.36;
    x_t = 0.36;
    rand_x_t = 0.93;
    alpha = 5.0;

    r_wOlambda_T_near = zeros(size(rand_r_omega));
    
    %genearted close wall r_w/lambda for rand_r_omega <= rand_x_t (Log-normal region)
    idx_log_normal = rand_r_omega <= rand_x_t;
    r_wOlambda_T_near(idx_log_normal) = exp(mu + sigma * norminv(rand_r_omega(idx_log_normal)));
    
    %genearted close wall r_w/lambda for rand_r_omega > rand_x_t (Power-law region)
    idx_power_law = rand_r_omega > rand_x_t;
    r_wOlambda_T_near(idx_power_law) = x_t * ((1 - rand_r_omega(idx_power_law)) / (1 - rand_x_t)).^(-1 / alpha);

    % Generate for far from the wall

    obj.ro_r_u_omega_far_wall = 0.45;
    n = 1e7; % Number of samples
    repro_numb = 143;
    [rand_u_omega, rand_r_omega] = genGaussCop(obj.ro_r_u_omega_far_wall,...
        n, repro_numb);

    % Generating u_omega
    a = 0.4;
    b = 4.5;
    m = 0.44;
    u_wOu_tau_far = -m*log(-rand_u_omega*(exp(-a/m)-exp(-b/m))+exp(-a/m));

    % Generating r_omega
    mu = -1.94;
    sigma = 0.36;
    x_t = 0.25;
    rand_x_t = 0.93;
    alpha = 4.5;

    r_wOlambda_T_far = zeros(size(rand_r_omega));
    
    %genearted close wall r_w/lambda For For rand_r_omega <= rand_x_t (Log-normal region)
    idx_log_normal = rand_r_omega <= rand_x_t;
    r_wOlambda_T_far(idx_log_normal) = exp(mu + sigma * norminv(rand_r_omega(idx_log_normal)));
    
    %genearted close wall r_w/lambda For For rand_r_omega > rand_x_t (Power-law region)
    idx_power_law = rand_r_omega > rand_x_t;
    r_wOlambda_T_far(idx_power_law) = x_t * ((1 - rand_r_omega(idx_power_law)) / (1 - rand_x_t)).^(-1 / alpha);

end