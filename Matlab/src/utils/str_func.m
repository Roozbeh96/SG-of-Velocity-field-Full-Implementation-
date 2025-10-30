function str_func(obj)
    %{
    More info about how to compute the structure function in:
    https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/
    article/stochastic-modal-velocity-field-in-roughwall-turbulence/
    800D97E60609F4FE4D92847CCDB1C7A0
    %}
    uprime = obj.Gen_u_HRVF - mean(obj.Gen_u_HRVF,3);
    concat_uprime = reshape(uprime,size(uprime, 1), []);
    win12del = ceil(12*obj.delta/(obj.Gen_x_HRVF(2)-obj.Gen_x_HRVF(1)));
    Numb_frame = floor(size(concat_uprime,2)/win12del);
    obj.D11=zeros(size(concat_uprime,1),win12del+1,Numb_frame);
    obj.r11=zeros(1,win12del+1);
    obj.epsilon_str = zeros(size(concat_uprime,1),Numb_frame);
    obj.eta_str = zeros(size(concat_uprime,1),Numb_frame);
    progressbar('Computing Structure Function')
    for r=0:win12del
        obj.r11(1,r+1)=(obj.Gen_x_HRVF(2)-obj.Gen_x_HRVF(1))*r;
    end
    
    for S = 1:Numb_frame
        for elev = 1:size(obj.z,1)
            for r=0:size(obj.r11,2)-1
                obj.D11(elev,r+1,S)=mean((concat_uprime(elev,(S-1)*win12del+1:S*win12del+1-r)-...
                    concat_uprime(elev,(S-1)*win12del+1+r:S*win12del+1)).^2,2);
            end
        end
        [max_val]=max(obj.D11(:,:,S).*obj.r11.^(-2/3),[],2);
        obj.epsilon_str(:,S) = (2./max_val).^(-3/2);
        obj.eta_str(:,S)=(obj.nu^3./obj.epsilon_str(:,S)).^0.25;
        progressbar((S)/Numb_frame)
    end

end