function [xnume, unume] = postprocessing(nele,coord,connect,un,xi);

% POST PROCESSING
% ---------------

xnume = []; u_dis = []; u_rot = [];
for el = 1:nele
    nd1 = connect(el,2); %node 1
    nd2 = connect(el,3); %node 2
    x_n = [coord(nd1,2), coord(nd2,2)];   %coordinate data of element
    u_n = un(2*nd1-1:2*nd2);    %Disp data at each node of the selected element
    Le = x_n(2) - x_n(1); %Element length
    Nx = [(1-xi)/2, (1+xi)/2];
    N1 = (2-3*xi+xi.^3)/4;
    N2 = (1-xi -xi.^2 +xi.^3)/4;
    N3 = (2 + 3*xi -xi.^3)/4;
    N4 = (-1 -xi + xi.^2 + xi.^3)/4;
    Nu = [N1, Le*N2/2, N3, Le*N4/2];
    
    xnume = [xnume;Nx*x_n'];
    u_dis = [u_dis;Nu*u_n];
    
    dN1 = (-3+3*xi.^2)/4;
    dN2 = (-1-2*xi+3*xi.^2)/4;
    dN3 = (3-3*xi.^2)/4;
    dN4 = (-1+2*xi+3*xi.^2)/4;
    
    Nu_rot = 2/Le*[dN1, Le*dN2/2, dN3, Le*dN4/2];
    u_rot = [u_rot;Nu_rot*u_n];
end
unume = [u_dis, u_rot];
