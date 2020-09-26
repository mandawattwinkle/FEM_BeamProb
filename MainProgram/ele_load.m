function F = ele_load(xi,qele,xvec)

Le = xvec(2) - xvec(1); %Element length
N = 1/4*[(1-xi)^2*(2+xi) ,Le/2*(1-xi)^2*(1+xi) , (1+xi)^2*(2-xi) , Le/2*(1+xi)^2*(xi-1)];
N1 = [(1-xi)/2 , (1+xi)/2];
x = N1*xvec';
q = qele(2) + qele(3)*x + qele(4)*x^2;  %Load distribution
F = N'*q*Le/2;
end