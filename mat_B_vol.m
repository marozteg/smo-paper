function B = mat_B_vol(young, elem, coord, dof, adof)
    ndof = max(max(dof));
    nel = size(elem,1);
    B = zeros(ndof,nel);
    cath = coord(elem(:,2),:) - coord(elem(:,1),:); % nel x 3
    length = vecnorm(cath'); % 1 x nel
    p = cath'./length; % (3 x nel)./(1 x nel) = 3 x nel
    coef = sqrt(young./(length'.^2)); % x=volume=AL, K(x) = sum(i=1, i=m, B(:,i)*B(:,i)'*x(i)); EA/L
    for i = 1:nel
        B(dof(i,:),i) = coef(i)*[p(:,i);-p(:,i)];
    end
    B = B(adof,:);
end