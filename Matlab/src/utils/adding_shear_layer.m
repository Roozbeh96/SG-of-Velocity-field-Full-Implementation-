function obj = adding_shear_layer(obj)
%{
 In this function, we smooth the sharp shear layer between the UMZs.
 detail are in section 2.2 of the paper:
 https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/
 stochastic-modal-velocity-field-in-roughwall-turbulence/
 800D97E60609F4FE4D92847CCDB1C7A0#article
%}

    delta_omega = 0.4*obj.lambda;
    theta = linspace(0, pi/2, 25);
    % One sample of step-like velocity profile. You can plot it before and
    % after smoothing shear-layer.
    %{
    figure
    plot(obj.Gen_u_prof(:,1),obj.z)
    %}
    progressbar('Adding Shear Layer')
    for i=1:size(obj.Gen_u_prof,2)
    
        Delta_u_m = diff(obj.Gen_u_prof(:,i));
        Delta_w_m = diff(obj.Gen_w_prof(:,i));
    
        change_indices = find(Delta_u_m ~= 0.0);
    
        for j=1:size(change_indices,1)
    
            Delta_u_mi=Delta_u_m(change_indices(j,1),1);
            Delta_w_mi=Delta_w_m(change_indices(j,1),1);
            
            u=obj.Gen_u_prof(change_indices(j,1),i)+Delta_u_mi/2*[(1-cos(theta)),(1+sin(theta))];
            w=obj.Gen_w_prof(change_indices(j,1),i)+Delta_w_mi/2*[(1-cos(theta)),(1+sin(theta))];
            z=(obj.z(1,change_indices(j,1))+obj.z(1,change_indices(j,1)+1))/2+...
                delta_omega/2*[(sin(theta)-1),(1-cos(theta))];
    
    
            [~,c]=find(z(1,1)<=obj.z & z(1,end)>=obj.z);
    
            for k=1:size(c,2)
    
                [~,cc]=find(obj.z(1,c(1,k))<=z,1,"first");
                
    
                obj.Gen_u_prof(c(1,k),i)=u(1,cc);
                obj.Gen_w_prof(c(1,k),i)=w(1,cc);
    
    
            end 
        end
        progressbar((i)/size(obj.Gen_u_prof,2))
    end
end