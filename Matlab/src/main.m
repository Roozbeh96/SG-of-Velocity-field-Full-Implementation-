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

figure
subplot(4,1,1)
plot(mean(Gen_sample.Gen_u_prof,2)/Gen_sample.u_tau,...
    Gen_sample.z/Gen_sample.z_0)
hold on
plot(1/kappa*log(1:1000),1:1000)
set(gca,'Yscale','log','TickLabelInterpreter','latex','FontSize',13)
xlabel('$U/u_{\tau}$','Interpreter','Latex','FontSize',14);
ylabel('$z/z_0$','Interpreter','Latex','FontSize',14);

subplot(4,1,2)
plot(Gen_sample.z/Gen_sample.z_0,...
    var(Gen_sample.Gen_u_prof,0,2)/Gen_sample.u_tau^2)
set(gca,'Xscale','log','TickLabelInterpreter','latex','FontSize',13)
xlabel('$z/z_0$','Interpreter','Latex','FontSize',14);
ylabel('$Var(u)/u_{\tau}^2$','Interpreter','Latex','FontSize',14);
ylim([0 7.5])

subplot(4,1,3)
plot(Gen_sample.z/Gen_sample.z_0,...
    var(Gen_sample.Gen_w_prof,0,2)/Gen_sample.u_tau^2)
set(gca,'Xscale','log','TickLabelInterpreter','latex','FontSize',13)
xlabel('$z/z_0$','Interpreter','Latex','FontSize',14);
ylabel('$Var(w)/u_{\tau}^2$','Interpreter','Latex','FontSize',14);
ylim([0 3])

subplot(4,1,4)
plot(Gen_sample.z/Gen_sample.z_0,...
    -sum((Gen_sample.Gen_u_prof - mean(Gen_sample.Gen_u_prof,2)) ...
    .*(Gen_sample.Gen_w_prof - mean(Gen_sample.Gen_w_prof,2)),2) ...
    / (Gen_sample.N_prof - 1)/Gen_sample.u_tau^2)
set(gca,'Xscale','log','TickLabelInterpreter','latex','FontSize',13)
xlabel('$z/z_0$','Interpreter','Latex','FontSize',14);
ylabel('$-Covar(u,w)/u_{\tau}^2$','Interpreter','Latex','FontSize',14);
ylim([0 1.5])
%% Stochastic Generation of Velocity Field (SGVF)

%{ 
    For the generated velocity field, you need to prescribe the distance
    between the velocity profiles which here is lambda. 
%}
Delta_x = lambda;
Gen_sample.SGVF(Delta_x)

%% Stochastic Generation of Vortex cores 

Gen_sample.SGVorX()




