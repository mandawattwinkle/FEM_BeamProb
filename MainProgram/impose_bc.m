function [K,F] = impose_bc(nele,K,F,BC_data)

% Modification of load vector 
% ---------------------------

supp_dof = [];

for i = 1:size(BC_data,1)
    nd = BC_data(i,1);
    dof = BC_data(i,2);
    value = BC_data(i,3);
    Gdof = 2*(nd-1)+dof;        %finding the global DOF of given BC
    supp_dof = [supp_dof,Gdof];
    if (value~=0)
        for j =1:2*(nele+1)
            F(j) = F(j) - K(j,Gdof)*value;  %imposition of non zero BC vale in load vector
        end
    end
end

supp_dof = sort(supp_dof);

% Reducing Global matrix (Imposition of Boundary condition)
% ---------------------------------------------------------

for i = 1:size(supp_dof,2)
    dof = supp_dof(i);  %Getting the Global DOF of Given BC
    if (dof == 1)       %Loop for reducing global matrix
        K = K(dof+1:end, dof+1:end);
        F = F(dof+1:end);
    elseif (dof == 2*(nele+1))
            K = K(1:dof-1, 1:dof-1);
            F = F(1:dof-1);
    else
        K = K([1:dof-1, dof+1:end],[1:dof-1, dof+1:end]);
        F = F([1:dof-1, dof+1:end]);
    end
    if (i~=size(supp_dof,2))    %for consiering the already reduced steps
        supp_dof(i+1:end) = supp_dof(i+1:end)-1;
    end
end
end