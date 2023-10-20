function dof = gs_dof(elem)
    % degree of freedom
    nel = size(elem,1);
    dof = zeros(nel,3*2);
    for i = 1:nel
        no1 = elem(i,1);
        no2 = elem(i,2);
        o1 = (no1-1)*3;
        o2 = (no2-1)*3;
        dof(i,:) = [o1+1 o1+2 o1+3 o2+1 o2+2 o2+3];
    end
end