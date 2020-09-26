function [K,F] = stiff_load(nele,ngauss,coord,connect,xivec,wvec,E,Ie,q_load);

%Elemental Stiffness matrix and Load vector
%------------------------------------------
K = zeros(2*(nele+1),2*(nele+1));   %Initialising Global stiffness and load matrixk
F = zeros(2*(nele+1),1);
for i = 1:nele  %Loop for each element calculation
    nd1 = connect(i,2); %node 1
    nd2 = connect(i,3); %node 2
    x = [coord(nd1,2), coord(nd2,2)];   %coordinate data of element
    Eele = E(i);    %Elementel youngs modulus
    Iele = Ie(i);   %Elemental moment of inertia
    if any(q_load(:,1)==i)  % checking wheather distributive load is given
        ii = find(q_load(:,1) ==i);
        qele = q_load(ii,:);    %getting the coefficients of distributive load for quadratic equation
    else
        qele = zeros(1,4);
    end
    
    
    Kele = zeros(4,4);
    Fele = zeros(4,1);
    for j =1:ngauss     %loop for guass quadrature
        xi = xivec(j);   %getting xi values
        w = wvec(j);     %corresponding w values
        Kele(1:4,1:4) = Kele(1:4,1:4) + ele_stiff(xi,Eele,Iele,x)*w;
        Fele(1:4) = Fele(1:4) + ele_load(xi,qele,x)*w;
    end
    
    % Assembly
    %---------
    vec  = [2*nd1-1, 2*nd1, 2*nd2-1, 2*nd2];    %Finding global DOF 

    for ii = 1:4
        for jj = 1:4
            K(vec(ii),vec(jj)) = K(vec(ii),vec(jj)) + Kele(ii,jj);  %Adding elemental stiffness matrix to Global stiffness matrix
        end
        F(vec(ii)) = F(vec(ii)) + Fele(ii); %%Adding elemental Load vvector to Global Load vector
    end
end
end