function F = loadcases(coord, load_cases, adof)

    nn = size(coord,1);
    nlc = size(load_cases,2);
    F = zeros(3*nn,nlc);
    for lc = 1:nlc
        loadmat = cell2mat(load_cases(lc));
        for i = 1:size(loadmat,1)
            for node = 1:nn
                if norm(coord(node,:) - loadmat(i,1:3))<1e-5
                    F((node-1)*3 + 1:(node-1)*3 + 3,lc) = loadmat(i,4:6)';
                end
            end
        end
    end
    F = F(adof,:);
end