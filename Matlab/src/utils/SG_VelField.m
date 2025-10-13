function [Gen_u_prof_reorg, Gen_w_prof_reorg, hist_corr] = SG_VelField(obj)
%SG_VELFIELD Summary of this function goes here
%   Detailed explanation goes here
    %If you do not want to smooth the shear layer in the generated
    %profiles, please comment this function out.
    %#ok<*AGROW>
    obj = adding_shear_layer(obj);
    %{
      Size of storage is Reynold's number dependent. 
    %}
    NP_storage=10^floor(log10(10*sqrt(obj.u_tau*obj.delta/obj.nu)));
    size_reorg_VField=size(obj.Gen_u_prof,2)-NP_storage;


    %{
      Since we want to sample W/O replacement, we take a copy of the 
        generated profiles, named it repository.  
    %}
    Gen_u_prof_repo = obj.Gen_u_prof;
    Gen_w_prof_repo = obj.Gen_w_prof;
    stepuprimeWT7 = obj.Gen_u_prof-mean(obj.Gen_u_prof,2);
    %{
      Sampled W/O replacement from repo, transferred to storage.  
    %}
    [u_storage,index]=datasample(Gen_u_prof_repo,NP_storage,2,'Replace',false);
    w_storage=Gen_w_prof_repo(:,index);
    uprimesampledatasto=stepuprimeWT7(:,index);
    
    Gen_u_prof_repo(:,index)=[];
    Gen_w_prof_repo(:,index)=[];
    stepuprimeWT7(:,index)=[];
    
    Gen_u_prof_reorg=zeros(size(obj.z,2),size_reorg_VField);
    Gen_w_prof_reorg=zeros(size(obj.z,2),size_reorg_VField);
    uprimefieldGenWT7=zeros(size(uprimeWT7Pro,1),size_reorg_VField);

    
    % sample first profile for the velocity field.
    r = randi([1, size(u_storage,2)]);
    Gen_u_prof_reorg(:,1)=u_storage(:,r);
    Gen_w_prof_reorg(:,1)=w_storage(:,r);
    uprimefieldGenWT7(:,1)=uprimesampledatasto(:,r);

    % Remove the selected profile from the storage
    u_storage(:,r)=[];
    w_storage(:,r)=[];
    uprimesampledatasto(:,r)=[];

    % Refill the storage with a randomly selected profile from repo.
    r = randi([1, size(Gen_u_prof_repo,2)]);
    u_storage=[u_storage,Gen_u_prof_repo(:,r)]; 
    w_storage=[w_storage,Gen_w_prof_repo(:,r)]; 
    uprimesampledatasto=[uprimesampledatasto,stepuprimeWT7(:,r)]; 

    % Remove the selected profile from the repo
    Gen_u_prof_repo(:,r)=[];
    Gen_w_prof_repo(:,r)=[];
    stepuprimeWT7(:,r)=[];

    C=zeros();
    hist_corr=cell(1,size(Gen_u_prof_reorg,2)-1);
    for prof_num_Vfield=2:size(Gen_u_prof_reorg,2)
        
        tic;
        for prof_num_storage=1:size(u_storage,2)
 
%             [c,~] = xcorr(ufieldGenWT7(:,j-1,k),usampledatasto(:,i),0,'normalized');
            [c,~] = xcorr(uprimefieldGenWT7(:,prof_num_Vfield-1),uprimesampledatasto(:,prof_num_storage),0,'normalized');
            C(1,prof_num_storage)=c;
            
        end
        [~,idx] = sort(C,'descend');
        C_sort=C(1,idx);

        % We can save the history of the correlation coefficient.
        
        hist_corr{1,prof_num_Vfield-1}=dictonary('cross-correlation with next profile',C_sort(1,1));
        
        % Filling the VField w/ the most correlated profile.
        Gen_u_prof_reorg(:,prof_num_Vfield)=u_storage(:,idx(1,1));
        Gen_w_prof_reorg(:,prof_num_Vfield)=w_storage(:,idx(1,1));
        uprimefieldGenWT7(:,prof_num_Vfield)=uprimesampledatasto(:,idx(1,1));

        % Remove the selected profile from the storage.
        u_storage(:,idx(1,1))=[];
        w_storage(:,idx(1,1))=[];
        uprimesampledatasto(:,idx(1,1))=[];
        
        % Refill the storage with a randomly selected profile from repo.
        r = randi([1, size(Gen_u_prof_repo,2)]);
        u_storage=[u_storage,Gen_u_prof_repo(:,r)]; 
        w_storage=[w_storage,Gen_w_prof_repo(:,r)]; 
        uprimesampledatasto=[uprimesampledatasto,stepuprimeWT7(:,r)]; 

        % Remove the selected profile from the repo
        Gen_u_prof_repo(:,r)=[];
        Gen_w_prof_repo(:,r)=[];
        stepuprimeWT7(:,r)=[];

        C=zeros();
        
        if mod(prof_num_Vfield, 500) == 0
            elapsedTime = toc;  % Elapsed time up to the checkpoint
            estimatedRemainingTime = (size(Gen_u_prof_reorg,2) - prof_num_Vfield) * (elapsedTime);
            fprintf('Estimated remaining time for reorganization: %.3f minutes\n', estimatedRemainingTime/60);
        end
    end
end

