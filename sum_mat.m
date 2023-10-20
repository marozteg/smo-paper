function M = sum_mat(nmin,nmax,M0,Mcel,x)
    M = M0;
    for i = nmin:nmax
        M = M + Mcel{i}*x(i);
    end
end