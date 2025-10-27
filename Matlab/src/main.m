%{
To Stochatically Generate Velocity Profiles (SGVP), reorganize, and add vortex
cores, you need to provide initial physical parameters such as
u_{tau}(friction velocity), z_{0}(aerodynamic roughness length),
delta(boundary layer thickness), lambda(Taylor micro-scale). The details 
about the rest of the parameter will be providing the following lines.
If you need details about the statistics, refer to the following papers:

**For SG of Velocity Profiles:
https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/
stochastic-modelling-of-the-instantaneous-velocity-profile-in-roughwall-
turbulent-boundary-layers/492F3CD03C8C3E7ED306E9117B848B5E

**For SG of Velocity field:
https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/
stochastic-modal-velocity-field-in-roughwall-turbulence/
800D97E60609F4FE4D92847CCDB1C7A0

**For SG of Vortex core:
coming soon
%}


root = fileparts(mfilename('fullpath'));  % project_root
addpath(fullfile(root, 'utils')); %addpath for using functions in utils.

%{
The physical properties of the wind tunnel (m1) are used. You can find 
details about this dataset in the papers.
%}
u_tau = 0.39; %[m/s]
z_0 = 6.2e-4; %[m]
kappa = 0.39;
delta = 0.4; %[m]
lambda = 0.01; %[m]
nu = 1.57e-5; %[m^2/s]
%{
Define the wall-normal extent of the generated profiles which covers
logarithmic layer. Log region starts from z+ = 30, continues till z/delta =
0.3. The values here are conservatively considered. The wall-normal
resolution of profiles is selected 0.4*lambda/10. Since average thickness
of shear layer is 0.4*lambda, I divided it into 10 to sesolve the shear
layer.
%}
z_min = 50*nu/u_tau; %[m]
z_max = 0.25*delta; %[m]
res_z = 0.4*lambda/10; %[m]
% Let us generate an object. 

Gen_sample = stochastic_generation(u_tau, z_0, delta, lambda, nu,z_min, z_max,...
    res_z);

%% Stochastic Generation of Velocity Profiles (SGVP)
%{
For the stochastic generation of velocity profiles, we need to start from
z_min, then stochastocally generate thickness (step), streamwise (or modal)
velocity, and vertical velocity. Repeat the process of step generation 
untill we reach to z_max, and then stop. We can independently generate 
profiles as much as we want. The minimum number of the profiles for having 
convergence in statistics (stationary dataset) is ~100 profiles. However, 
since we want to generate long velocity field in the next step, we need 
more profiles.
The parameter we need here is the Pearson correlation coefficient between
streamwise velocity and vertical velocity components. This parameter is
necessary in Reynolds shear stress (covariance) reproduction. This value in
the wind tunnel (m1) dataset is -0.4. You can tweak it w.r.t. your dataset.
N_prof is number of profiles we want to generate.
%}

ro_uw = -0.4;
N_prof = 1e4;
Gen_sample.SGVP(ro_uw, N_prof)

%% Stochastic Generation of Velocity Field (SGVF)

%{ 
    For the generated velocity field, you need to prescribe the distance
    between the velocity profiles which here is lambda. 
%}
Delta_x = lambda;
Gen_sample.SGVF(Delta_x)

%% Increasing the resolution of the SGVF and compute the lambda_ci field
%{
 At the first step, we need to divide the VF into 1 * delta (L = 1)
 intervals, and increase the resolution. For increasing the resolution,
 the factor that should be increased is another input parameter. Here 
 Kr is the factor parameter which is consistent with the parameter in 
 the paper.  All these steps are done increasing_res function.
%}
L = 1; 
Kr = 32;
increasing_res(Gen_sample, L, Kr);
% In the next step, we need to find the Lambda_ci field to add the
% vortices.
Lambda_ci(Gen_sample)

%% Generating radius, maximum azimuthal velocity, and corrcoeff(u,w)
%{
    At this step, we genrate radius(r_{\omega}), maximum azimuthal velocity
    (u_{\omega}), and the correlation coefficient between u and w component
    of the generated vortex \rho^{VorX}_{u,w}.
%}

delUwOu_tau_near = load('delUwOu_tau_near.mat').delUwOu_tau_near;
delUwOu_tau_Far = load('delUwOu_tau_Far.mat').delUwOu_tau_Far;
rwOlambda_T_near = load('rwOlambda_T_near.mat').rwOlambda_T_near;
rwOlambda_T_Far = load('rwOlambda_T_Far.mat').rwOlambda_T_Far;
%% Stochastic Generation of Vortex cores
Gen_sample.SGVorX()




