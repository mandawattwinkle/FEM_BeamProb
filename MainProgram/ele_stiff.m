function K = ele_stiff(xi,E,I,x)
Le = x(2) - x(1);    %element Length
B = 4*[3*xi/2, Le/4*((3*xi)-1), -3*xi/2, Le/4*((3*xi)+1)]/Le^2;

K = E*I*Le*B'*B/2;  %element stiffness
end