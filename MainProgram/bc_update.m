function un = bc_update(ureduce,BC_data)

% Full Solution vector (Free + Prescribed D. O. F.)
% -------------------------------------------------

un = ureduce;
for i = 1:size(BC_data,1)
    nd = BC_data(i,1);
    dof = BC_data(i,2);
    value = BC_data(i,3);   %Getting Value of the prescribed DOF
    Gdof = 2*(nd-1)+dof;    %Finding global DOF
    if (Gdof ==1)           %Loop for adding prescribed DOF in the Global Disp Vector
        un = [value;un];
    else
        un = [un(1:Gdof-1);value;un(Gdof:end)];
    end
end