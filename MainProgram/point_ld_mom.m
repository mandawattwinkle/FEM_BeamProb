function F = point_ld_mom(F,P_load,P_moment)
% Point Load and Point Moment
% ---------------------------

for i = 1:size(P_load,1)
    nd = P_load(i,1);    %getting node in which point load acts
    F_point = P_load(i,2);  %Point load acts on the node
    F(2*nd-1) = F(2*nd-1)+F_point;
end

for i = 1:size(P_moment,1)
    nd = P_moment(i,1);    %getting node in which point Moment acts
    M_point = P_moment(i,2);  %Point load acts on the node
    F(2*nd) = F(2*nd)+M_point;
end
end
