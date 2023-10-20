function Mdx = mat_Mdx_vol(rho,elem,coord,dof,adof)
    % derivative of M relative to x=volume 
    nel = size(elem,1);
    Mdx = cell(nel,1);
    ndof = max(max(dof));
    r = coord(elem(:,2),:) - coord(elem(:,1),:);
    length = vecnorm(r')';
    coef = rho.*length/6;
    for i = 1:nel
        MatElem = coef(i)*( diag(2*ones(6,1)) + diag(ones(3,1),3) + diag(ones(3,1),-3) );
        Matrix = zeros(ndof,ndof);
        gl12 = dof(i,:);
        Matrix(gl12,gl12) = MatElem;
        Mdx{i} = sparse(Matrix(adof,adof));
    end
end
