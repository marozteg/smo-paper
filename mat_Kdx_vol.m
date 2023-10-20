function Kdx = mat_Kdx_vol(young,elem,coord,dof,adof)
    B = mat_B_vol(young, elem, coord, dof, adof);
    nel = size(elem,1);
    Kdx = cell(nel,1);
    for i = 1:nel
        Kdx{i} = B(:,i)*B(:,i)';
    end
end