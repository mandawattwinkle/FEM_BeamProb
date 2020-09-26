function [x_ana,u,theeta] = analytical(E0,Ie0,q0,F0,L1,L2)


x_ana=[0:.05:L1+L2];
theeta =zeros(1,length(x_ana));
u = zeros(1,length(x_ana));
for i = 1:length(x_ana)
    x=x_ana(i);
    if (x<L1)
        c1 = +(F0*L1^2/2)+(q0*L1^3/60);
        c2 = -(F0*L1^3/6)-(q0*L1^4/360);
        theeta(i) = 1/(E0*Ie0)*((F0*((L2*x)-((x-L1)^2/2)))+(q0*(x-L1)^5/(60*L1^2))+c1);
        u(i) = 1/(E0*Ie0)*((F0*((L2*x^2/2)-((x-L1)^3/6)))+(q0*(x-L1)^6/(L1^2*360))+(c1*x)+c2);
    else
        c1 = (q0*L1^3/60);
        c2 = -(q0*L1^4/360);
        theeta(i) = 1/(E0*Ie0)*((F0*((L2*x)+(L1*x)-(x^2/2)))+c1);
        u(i) = 1/(E0*Ie0)*((F0*((L2*x^2/2)+(L1*x^2/2)-(x^3/6)))+(c1*x)+c2);
    end
end







