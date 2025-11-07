%% Load data
load("Gen_sample.mat")

%% plot first and second moments of statistics Generated profiles
kappa = 0.39;
figure
subplot(4,1,1)
plot(mean(Gen_sample.Gen_u_prof,2)/Gen_sample.u_tau,...
    Gen_sample.z/Gen_sample.z_0,'Color',[0,0,1],'linewidth',2)
hold on
plot(1/kappa*log(1:1000),1:1000,'Color',[1,0,1],'linewidth',1.5)
set(gca,'Yscale','log','TickLabelInterpreter','latex','FontSize',13)
xlabel('$\mathrm{U}/u_{\tau}$','Interpreter','Latex','FontSize',14);
ylabel('$z/z_0$','Interpreter','Latex','FontSize',14);

subplot(4,1,2)
plot(Gen_sample.z/Gen_sample.z_0,...
    var(Gen_sample.Gen_u_prof,0,2)/Gen_sample.u_tau^2,'Color',[0,0,1],...
    'linewidth',2)
set(gca,'Xscale','log','TickLabelInterpreter','latex','FontSize',13)
xlabel('$z/z_0$','Interpreter','Latex','FontSize',14);
ylabel('$\mathrm{Var(u)}/u_{\tau}^2$','Interpreter','Latex','FontSize',14);
ylim([0 7.5])

subplot(4,1,3)
plot(Gen_sample.z/Gen_sample.z_0,...
    var(Gen_sample.Gen_w_prof,0,2)/Gen_sample.u_tau^2,'Color',[0,0,1],...
    'linewidth',2)
set(gca,'Xscale','log','TickLabelInterpreter','latex','FontSize',13)
xlabel('$z/z_0$','Interpreter','Latex','FontSize',14);
ylabel('$\mathrm{Var(w)}/u_{\tau}^2$','Interpreter','Latex','FontSize',14);
ylim([0 3])

subplot(4,1,4)
plot(Gen_sample.z/Gen_sample.z_0,...
    -sum((Gen_sample.Gen_u_prof - mean(Gen_sample.Gen_u_prof,2)) ...
    .*(Gen_sample.Gen_w_prof - mean(Gen_sample.Gen_w_prof,2)),2) ...
    / (Gen_sample.N_prof - 1)/Gen_sample.u_tau^2,'Color',[0,0,1],...
    'linewidth',2)
set(gca,'Xscale','log','TickLabelInterpreter','latex','FontSize',13)
xlabel('$z/z_0$','Interpreter','Latex','FontSize',14);
ylabel('$-\mathrm{Covar(u,w)}/u_{\tau}^2$','Interpreter','Latex','FontSize',14);
ylim([0 1.5])


%% plot first and second moments of statistics Generated Velocity field
figure
subplot(4,1,1)
plot(mean(mean(Gen_sample.Gen_u_HRVF,3),2)/Gen_sample.u_tau,...
    Gen_sample.z/Gen_sample.z_0,'Color',[1,0,0],'linewidth',2)
hold on
plot(1/kappa*log(1:1000),1:1000,'Color',[1,0,1],'linewidth',1.5)
set(gca,'Yscale','log','TickLabelInterpreter','latex','FontSize',13)
xlabel('$\mathrm{U}/u_{\tau}$','Interpreter','Latex','FontSize',14);
ylabel('$z/z_0$','Interpreter','Latex','FontSize',14);

subplot(4,1,2)
plot(Gen_sample.z/Gen_sample.z_0,...
    mean(var(Gen_sample.Gen_u_HRVF,0,3),2)/Gen_sample.u_tau^2,...
    'Color',[1,0,0],'linewidth',2)
set(gca,'Xscale','log','TickLabelInterpreter','latex','FontSize',13)
xlabel('$z/z_0$','Interpreter','Latex','FontSize',14);
ylabel('$\mathrm{Var(u)}/u_{\tau}^2$','Interpreter','Latex','FontSize',14);
ylim([0 7.5])

subplot(4,1,3)
plot(Gen_sample.z/Gen_sample.z_0,...
    mean(var(Gen_sample.Gen_w_HRVF,0,3),2)/Gen_sample.u_tau^2,...
    'Color',[1,0,0],'linewidth',2)
set(gca,'Xscale','log','TickLabelInterpreter','latex','FontSize',13)
xlabel('$z/z_0$','Interpreter','Latex','FontSize',14);
ylabel('$\mathrm{Var(w)}/u_{\tau}^2$','Interpreter','Latex','FontSize',14);
ylim([0 3])

subplot(4,1,4)
plot(Gen_sample.z/Gen_sample.z_0,...
    mean(-sum((Gen_sample.Gen_u_HRVF - mean(Gen_sample.Gen_u_HRVF,3)) ...
    .*(Gen_sample.Gen_w_HRVF - mean(Gen_sample.Gen_w_HRVF,3)),3) ...
    / (size(Gen_sample.Gen_u_HRVF,3) - 1),2)/Gen_sample.u_tau^2,...
    'Color',[1,0,0],'linewidth',2)
set(gca,'Xscale','log','TickLabelInterpreter','latex','FontSize',13)
xlabel('$z/z_0$','Interpreter','Latex','FontSize',14);
ylabel('$-\mathrm{Covar(u,w)}/u_{\tau}^2$','Interpreter','Latex','FontSize',14);
ylim([0 1.5])
%% Spectral analysis

ind_ = 20;
ind__ = 96;
ind___ = 196;

figure
set(gcf,'Position',[928,562,407,383])
axes('Position',[0.181818181818182,0.143603133159269,0.778869778869779,0.827676240208877])
loglog((3e-4:1e-5:1),(3e-4:1e-5:1).^(-5/3),'linewidth',2,'color',[1.00,0.74,0.17],'lineStyle',':')
hold on
loglog(linspace(5e-7,5,size(Gen_sample.wavenumb,1)),ones(size(Gen_sample.wavenumb)),'k--','linewidth',1.5)
loglog(ones(size(Gen_sample.wavenumb)),linspace(1e-8,10,size(Gen_sample.wavenumb,1)),'k--','linewidth',1.5)


numPoints = 30; % Adjust as needed
logIndices = round(logspace(0, log10(size(Gen_sample.wavenumb, 1)), numPoints));

loglog(Gen_sample.wavenumb(logIndices,1)*mean(Gen_sample.eta_str(ind_,:),2),Gen_sample.PowSpecDen_k(ind_,logIndices)/...
    (mean(Gen_sample.epsilon_str(ind_,:),2)*Gen_sample.nu^5)^(1/4),'-.',...
    'linewidth',2,'color',[1.0, 0.0, 0.0])

loglog(Gen_sample.wavenumb(logIndices,1)*mean(Gen_sample.eta_str(ind__,:),2),Gen_sample.PowSpecDen_k(ind__,logIndices)/...
    (mean(Gen_sample.epsilon_str(ind__,:),2)*Gen_sample.nu^5)^(1/4),'-',...
    'linewidth',2,'color',[0.0, 1.0, 0.0])

loglog(Gen_sample.wavenumb(logIndices,1)*mean(Gen_sample.eta_str(ind___,:),2),Gen_sample.PowSpecDen_k(ind___,logIndices)/...
    (mean(Gen_sample.epsilon_str(ind___,:),2)*Gen_sample.nu^5)^(1/4),'--',...
    'linewidth',2,'color',[0.0, 0.0, 1.0])

xlabel('$\mathrm{k}_{1}\eta$','Interpreter','latex')
ylabel('$\mathrm{E}_{11}(\mathrm{k}_{1})/(\varepsilon \nu^5)^{1/4}$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',13,'XGrid','on','YGrid','on',...
    'YTick', [1e-6,1e0,1e6])
xlim([1e-4 10])
ylim([1e-6 1e6])
annotation(gcf,'textbox',...
    [0.328276792088245,0.885117493472585,0.093078536700127,0.082433208914058],...
    'String',{'$(\mathrm{k}_{1}\eta)^{-5/3}$'},...
    'Interpreter','latex',...
    'FontSize',13,...
    'FitBoxToText','off',...
    'EdgeColor','none');
legend('','','',...
    sprintf('HR+VorX-Gen(m1) z$/\\delta$ = %.2f',Gen_sample.z(ind_)/Gen_sample.delta),...
    sprintf('HR+VorX-Gen(m1) z$/\\delta$ = %.2f',Gen_sample.z(ind__)/Gen_sample.delta),...
    sprintf('HR+VorX-Gen(m1) z$/\\delta$ = %.2f',Gen_sample.z(ind___)/Gen_sample.delta),...
    'Interpreter','latex','FontSize',10,'Position',...
    [0.199385032377767,0.210296631208215,0.535213977832512,0.129128956259148],...
    'Numcolumns',1,'Orientation','vertical','color','none');
set(gca,'TickLabelInterpreter','latex','FontSize',13)
axis square