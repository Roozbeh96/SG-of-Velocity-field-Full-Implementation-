function MVG_kernel(obj, kern_xsize, kern_ysize, Sigma)

    mu = [0 0];

    [X,Y] = meshgrid(linspace(-kern_xsize/2 * obj.lambda,...
        kern_xsize/2 * obj.lambda,...
        ceil(kern_xsize*obj.lambda/(obj.Gen_x_HRVF(2)-obj.Gen_x_HRVF(1)))+1),...
        linspace(-kern_ysize/2 * obj.lambda,...
        kern_ysize/2 * obj.lambda,...
        ceil(kern_ysize*obj.lambda/(obj.z(2)-obj.z(1)))+1));
    kernel = mvnpdf([X(:) Y(:)],mu,Sigma);
    kernel = reshape(kernel,size(X));
    
    figure
    set(gcf,'Position',[806,788,560,271])
    contourf(X/obj.lambda, Y/obj.lambda, kernel, 100,'LineStyle','none')
    colorbar
    hcb2=colorbar;
    hcb2.TickLabelInterpreter = 'latex';
    set(gca,'TickLabelInterpreter','latex','FontSize',13)
    xlabel('x/$\lambda_{T}$','Interpreter','Latex','FontSize',14);
    ylabel('z/$\lambda_{T}$','Interpreter','Latex','FontSize',14);
    axis equal

    normalization_matrix = conv2(ones(size(obj.Gen_u_HRVF,1),size(obj.Gen_u_HRVF,2)),...
        kernel, 'same');
    for S =1:size(obj.HRVFu,3)
    
        obj.Gen_u_HRVF(:,:,S) = conv2(obj.Gen_u_HRVF(:,:,S), kernel, 'same')./normalization_matrix;
        obj.Gen_w_HRVF(:,:,S) = conv2(obj.Gen_w_HRVF(:,:,S), kernel, 'same')./normalization_matrix;

    end

end