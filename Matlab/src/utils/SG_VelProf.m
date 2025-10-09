function Gen_Prof = SG_VelProf(obj)
%SGVF generate step like velocity profiles. Detail of line will explain in
%the folllowing lines.

    [rand_u,rand_w] = genGaussCop(obj.ro_uw, 1e6);
    Gen_Prof = zeros(size(obj.z,2),obj.N_prof);
    
    
    for ind = 1:obj.N_prof
        
    
    end


end
