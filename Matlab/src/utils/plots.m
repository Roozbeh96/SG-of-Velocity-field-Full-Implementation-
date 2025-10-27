%% plot first and second moments of statistics
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

%%