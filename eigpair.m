function [delta, v] = eigpair(A,B,sense)
% given A and B positive semidefinite and sense
% compute:
%   if sense == "min"
%      delta = argmax{delta : B - delta*A positive semidefinite} 
%            = argmin{delta : Bv = delta*Av, v notin kerA}
%   if sense == "max"
%      delta = argmin{delta : delta*A - B positive semidefinite}
%            = argmax{delta : Bv = delta*Av, v notin kerA}
%   if A is positive definite,
%       v, the eigenvector associated to delta

    P = orth(full(A));
    lambda = eig(P'*B*P,P'*A*P);
    if nargout >= 1
        if sense == "max"
            delta = max(lambda);
        elseif sense == "min"
            delta = min(lambda);
        end
    end
    if nargout == 2
        if min(eig(A)) > 1e-7
            [V,D] = eig(B,A);
            [~,idx] = sort(diag(D));
            if sense == "max"
                v = V(:,idx(end));
            elseif sense == "min"
                v = V(:,idx(1));
            end
        end
    end
end