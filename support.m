function adof = support(coord, support_coord, kind)

    nn = size(coord,1);
    nsc = size(support_coord);
    fdof = [];
    for i = 1:nsc
        for node = 1:nn
            if norm(coord(node,:) - support_coord(i,:))<1e-5
                if kind(i) == "x"
                    fdof = [fdof (node-1)*3+1];
                elseif kind(i) == "y"
                    fdof = [fdof (node-1)*3+2];
                elseif kind(i) == "z"
                    fdof = [fdof (node-1)*3+3];
                elseif kind(i) == "xy"
                    fdof = [fdof (node-1)*3+1 (node-1)*3+2];
                elseif kind(i) == "xz"
                    fdof = [fdof (node-1)*3+1 (node-1)*3+3];
                elseif kind(i) == "yz"
                    fdof = [fdof (node-1)*3+2 (node-1)*3+3];
                elseif kind(i) == "xyz"
                    fdof = [fdof (node-1)*3+1 (node-1)*3+2 (node-1)*3+3];
                end
            end
        end
    end
    fdof = unique(fdof);
    adof = 1:3*nn;
    adof(fdof) = [];
end