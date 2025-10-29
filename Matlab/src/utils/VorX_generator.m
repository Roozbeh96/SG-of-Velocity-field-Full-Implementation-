function [u_wOu_tau_near, r_wOlambda_T_near, rho_uw_VorX_near,...
    u_wOu_tau_far, r_wOlambda_T_far, rho_uw_VorX_far]= VorX_generator(obj)

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

    % Check to plot the p.d.f. You can try it for other parameters.
    % histogram(r_wOlambda_T_near,'Normalization','pdf','NumBins',500)
    % [f, xi] = ksdensity(r_wOlambda_T_near);
    % plot(xi, f, 'r', 'LineWidth', 2)

    % Generating \rho^{VorX}_{u,w}
    C1 = 1.23;
    omega = 0.58;
    kesi = -0.61;
    beta = 2.37;

    pdf_skewedGauss = @(x) ...
    C1 * (1 - x.^2) .* (2 ./ omega .* normpdf((x - kesi) ./ omega) .* normcdf(beta .* (x - kesi) ./ omega));

    x_vals = linspace(-1, 1, 1000); 
    % Step 1: Calculate the p.d.f. over the range
    pdf_vals = pdf_skewedGauss(x_vals);
    
    % Step 2: Compute the c.d.f. using cumulative trapezoidal integration
    cdf_vals = cumtrapz(x_vals, pdf_vals);
    cdf_vals = cdf_vals ./ max(cdf_vals);  % Normalize to ensure c.d.f. ends at 1
    
    % Step 3: Invert the c.d.f. using interpolation
    inversecdf = @(rand_) interp1(cdf_vals, x_vals, rand_, 'linear', 'extrap');
    
    % Step 4: Generate uniform random numbers
    rng(18)
    rand_rho_uw_VorX = rand(n, 1);  % Uniform random numbers in [0, 1]
    
    % Step 5: Generate samples using the inverse c.d.f.
    rho_uw_VorX_near = inversecdf(rand_rho_uw_VorX);


    % Generate vortex parameters for far from the wall

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

    % Generating \rho^{VorX}_{u,w}
    C1 = 1.28;
    omega = 0.66;
    kesi = -0.54;
    beta = 1.80;

    pdf_skewedGauss = @(x) ...
    C1 * (1 - x.^2) .* (2 ./ omega .* normpdf((x - kesi) ./ omega) .* normcdf(beta .* (x - kesi) ./ omega));

    x_vals = linspace(-1, 1, 1000); 
    pdf_vals = pdf_skewedGauss(x_vals);
    cdf_vals = cumtrapz(x_vals, pdf_vals);
    cdf_vals = cdf_vals ./ max(cdf_vals);  
    inversecdf = @(rand_) interp1(cdf_vals, x_vals, rand_, 'linear', 'extrap');
    rng(61)
    rand_rho_uw_VorX = rand(n, 1);  
    rho_uw_VorX_far = inversecdf(rand_rho_uw_VorX);

end