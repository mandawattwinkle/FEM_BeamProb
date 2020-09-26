% opening a File named Result for printing the result of the Program
% ==================================================================

fid = fopen('../example c/Result.txt','w');
fprintf(fid,'Beam Problem\n');

fclose(fid);
% ===================
% Reading Input File:
% ==================
for ii = 1:2

str1 = sprintf(['../example c/input_file_beam_prob_c_',int2str(ii)]);

run(str1)


% Global Stiffnes Matrix and Global load vector
% ---------------------------------------------
% Function "stiff_load" calculate Global Stiffness Matrix and 
%     Global load vector due to distributed load
% ----------
% I N P U T
% = = = = = 
% nele   = No. of elements 
% ngauss = No. of gauss points for integration
% coord  = Nodal coordinates    % First Column is Node numbers
%                                 Second Column is Co-ordinate
% connect = Nodal Connectivities    % First Column is element number    
                                    % Second & Third Column are Nodes (in sequence)  
                                    % For that element. 
% xivec = Gauss points
% wvec  = weights
% E = Young's Modulus of the element
% Ie = Area Moment of inertia of the element
% q0 = Maximum distributed load (Triangular) magnitude
% L  = Length
%
% ----------
% O U T P U T
% = = = = =
% K = Global stiffness matrix
% F = Global load vector
%

[K,F] = stiff_load(nele,ngauss,coord,connect,xivec,wvec,E,Ie,q_load);

% Point load and Point moment
% ---------------------------
% This function "point_ld_mom" update Global load vector after incorporating point load 
% and point moment data
%
% ----------
% I N P U T
% = = = = = 
% F      = Global load vector before implementing point load and point moment
% P_load = Point load data  % First Column is Node number
                            % Second Column is Point load value
% P_moment = Point moment data       % First Column is Node number
                                     % Second Column is Point moment value
% -------------
% O U T P U T
% = = = = =====
% F = Global load vector after implementing point load and point moment

F = point_ld_mom(F,P_load,P_moment);


% ===========================
% Imposition of B.C.
K_glob = K;
F_glob = F;

% This function "impose_bc" update Global stiffness matrix and Global load vector 
% after incorporating boundary condition data
%
% ----------
% I N P U T
% = = = = = 
% K       = Global stiffness matrix before implementing Boundary condition data
% F       = Global load vector before implementing Boundary condition data
% BC_data = Boundary condition data        % First Column is Node number
                                           % Second Column is the prescribed D.O.F
                                           % Third Column is value of the prescribed D.O.F 
%
% -------------
% O U T P U T
% = = = = =====
% K       = Global stiffness matrix after implementing Boundary condition data
% F       = Global load vector after implementing Boundary condition data

[K,F] = impose_bc(nele,K,F,BC_data);


% Finding Solution
ureduce = inv(K)*F;

% Full Solution vector (Free + Prescribed D. O. F.)
% -------------------------------------------------
% This function "bc_update" update solution vector with values of prescribed DOFs
% 
% I N P U T
% =========
% ureduce = Solution vector just after inversion
%           It contains only free DOFs
% BC_data = Boundary condition data        % First Column is Node number
                                           % Second Column is the prescribed D.O.F
                                           % Third Column is value of the prescribed D.O.F 
% -------------
% O U T P U T
% = = = = =====
% un = Full solution vectors with Free and Prescribed DOF values
un = bc_update(ureduce,BC_data);

% Finding Reaction Force
Freac = K_glob*un;

% Post Processing: FEM displacement
xi = [-1:0.2:1]';          % Distribution of data points

% This function "postprocessing" calculate variable u at diffent distributed points across
% element from nodal values of u
% This function calculate variable u at diffent distributed points across
% element from nodal values of u
% ----------
% I N P U T
% = = = = = 
% nele   = No. of elements
% coord  = Nodal coordinates    % First Column is Node numbers
%                                 Second Column is Co-ordinate
% connect = Nodal Connectivities    % First Column is element number    
                                    % Second & Third Column are Nodes (in sequence)  
                                    % For that element. 
% xi = Points distributed for an element in master domain
% un = Nodal values of u
%
% ----------
% O U T P U T
% = = = = =
% xnume = x coordinates of the distributed points
% unume = values of u at distributed points

[xnume, unume] = postprocessing(nele,coord,connect,un,xi);

% PLOTTING
% ========

h = figure(1);

subplot(1,2,1)
set(gcf, 'Position', get(0,'Screensize'));
plot(xnume,unume(:,1));
grid on;
xlabel (' x (m)','fontsize',18);
ylabel (' w (m)','fontsize',18);
title ('Deflection curve','fontsize',20);   
hold on
subplot(1,2,2)
plot(xnume, unume(:,2));
grid on;
xlabel (' x (m)','fontsize',18);
ylabel (' \theta (rad)','fontsize',18);
title ('Rotation curve','fontsize',20);
hold on


% Printing results to output file
% ===============================

fid = fopen('../example c/Result.txt','a');

fprintf(fid,'\n\n===========================================================\n');
fprintf(fid,'E = %12.4e\nI1 = %12.4e;\tI2 = %12.4e\n\n',E0,I1,I2);
fprintf(fid,'Number of Elements used = %d\n',nele);
fprintf(fid,'Global Stiffness Matrix ''K'' is\n\n');
for i = 1:2*(nele+1)
    for j = 1:2*(nele+1)
    fprintf(fid,'%2.3e\t',K_glob(i,j));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n\nGlobal Load Vector ''F'' is\n\n');
for i = 1:2*(nele+1)
    fprintf(fid,'%14.4e\n',F_glob(i));
end
fprintf(fid,'\nReduced Stiffness Matrix after imposing Boundary Condition is\n\n');
for i = 1:size(K,1)
    for j = 1:size(K,1)
    fprintf(fid,'%2.3e\t',K(i,j));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\nReduced Load vector after imposing Boundary Condition is\n\n');
for i = 1:size(F)
    fprintf(fid,'%14.4e\n',F(i));
end
fprintf(fid,'\n Displacement vector\n');
for i = 1:size(un,1)
    fprintf(fid,'%14.4e\n',un(i));
end

fprintf(fid,'\nNode\t\tDisplacement\t\tRotation\n');
for i = 1:nele+1
    fprintf(fid,'%d\t\t\t%f\t\t\t%f\n',i,un(2*i-1),un(2*i));
end
    
fclose(fid);

end
subplot(1,2,1);
legend('No.of Elements = 3','No of Elements = 6');
hold on
subplot(1,2,2);
legend('No.of Elements = 3','No of Elements = 6');
hold on
saveas(h,'../example c/output_graph','fig');