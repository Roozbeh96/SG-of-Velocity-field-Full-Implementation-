function obj = SG_VelField(obj)
%SG_VELFIELD Summary of this function goes here
%   Detailed explanation goes here
    %If you do not want to smooth the shear layer in the generated
    %profiles, please comment this function out.
    obj = adding_shear_layer(obj);

    NPpool=1000;
    NP=61000;
    Nf=1;
    checkpoint_Interval = 500; 
    %{
        
    %}
    stepuWT7 = obj.Gen_u_prof;
    stepwWT7 = obj.Gen_w_prof;
    stepuprimeWT7 = obj.Gen_u_prof-mean(obj.Gen_u_prof,2);
    
    [usampledatasto,index]=datasample(stepuWT7,NPpool,2,'Replace',false);
    wsampledatasto=stepwWT7(:,index);
    uprimesampledatasto=stepuprimeWT7(:,index);
    
    stepuWT7(:,index)=[];
    stepwWT7(:,index)=[];
    stepuprimeWT7(:,index)=[];
    
    ufieldGenWT7=zeros(size(meanuWT7Pro,1),NP,Nf);
    wfieldGenWT7=zeros(size(meanwWT7Pro,1),NP,Nf);
    uprimefieldGenWT7=zeros(size(uprimeWT7Pro,1),NP,Nf);
    
    
       
        
    for k=1:Nf
        
    
        C=zeros();
        
        r = randi([1, size(usampledatasto,2)]);
        
        ufieldGenWT7(:,1,k)=usampledatasto(:,r);
        wfieldGenWT7(:,1,k)=wsampledatasto(:,r);
        uprimefieldGenWT7(:,1,k)=uprimesampledatasto(:,r);
        usampledatasto(:,r)=[];
        wsampledatasto(:,r)=[];
        uprimesampledatasto(:,r)=[];
        
        r = randi([1, size(stepuWT7,2)]);
        
        usampledatasto=[usampledatasto,stepuWT7(:,r)]; %#ok<AGROW>
        wsampledatasto=[wsampledatasto,stepwWT7(:,r)]; %#ok<AGROW>
        uprimesampledatasto=[uprimesampledatasto,stepuprimeWT7(:,r)]; %#ok<AGROW>
    
        stepuWT7(:,r)=[];
        stepwWT7(:,r)=[];
        stepuprimeWT7(:,r)=[];
        
        CorrelationGenWT7=zeros(1,size(ufieldGenWT7,2)-1);
        for j=2:size(ufieldGenWT7,2)
            
            tic;
            for i=1:size(usampledatasto,2)
     
    %             [c,~] = xcorr(ufieldGenWT7(:,j-1,k),usampledatasto(:,i),0,'normalized');
                [c,~] = xcorr(uprimefieldGenWT7(:,j-1,k),uprimesampledatasto(:,i),0,'normalized');
                C(1,i)=c;
                
    %                         C(1,i)=sum(uprimefieldGenWT7(:,j-1,k).*uprimesampledatasto(:,i),1)/...
    %                             ((size(uprimefieldGenWT7,1))*rms(uprimefieldGenWT7(:,j-1,k))*rms(uprimesampledatasto(:,i)));
                
            end
            [~,idx] = sort(C,'descend');
            Csort=C(1,idx);
            CorrelationGenWT7(1,j-1)=Csort(1,1);
            ufieldGenWT7(:,j,k)=usampledatasto(:,idx(1,1));
            wfieldGenWT7(:,j,k)=wsampledatasto(:,idx(1,1));
            uprimefieldGenWT7(:,j,k)=uprimesampledatasto(:,idx(1,1));%
            usampledatasto(:,idx(1,1))=[];
            wsampledatasto(:,idx(1,1))=[];
            uprimesampledatasto(:,idx(1,1))=[];
            
            r = randi([1, size(stepuWT7,2)]);
            usampledatasto=[usampledatasto,stepuWT7(:,r)]; %#ok<AGROW>
            wsampledatasto=[wsampledatasto,stepwWT7(:,r)]; %#ok<AGROW>
            uprimesampledatasto=[uprimesampledatasto,stepuprimeWT7(:,r)]; %#ok<AGROW>
    
            stepuWT7(:,r)=[];
            stepwWT7(:,r)=[];
            stepuprimeWT7(:,r)=[];
            C=zeros();
            
            
            
            if mod(j, checkpoint_Interval) == 0
                elapsedTime = toc;  % Elapsed time up to the checkpoint
                estimatedRemainingTime = (size(ufieldGenWT7,2) - j) * (elapsedTime);
                fprintf('Estimated remaining time (WT(m1)): %.3f minutes\n', estimatedRemainingTime/60);
            end
        end
    end
end

