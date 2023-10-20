function draw_nonstructmass(coord,adofmass,adof,y,minL)
    nn = size(coord,1);
    ng = size(y,1);
    for g = 1:ng
        coef = (minL/10)*2 / ((max(y{g}))^(1/3));
        nmg = size(adofmass{g},1);
        for i = 1:nn
            for dof = (i-1)*3+1:(i-1)*3+3
                for k = 1:nmg
                    adofmass_mat = adof(adofmass{g}{k}==1);
                    if ~isempty(find(adofmass_mat == dof, 1))
                        mmm = y{g}(k);
                        if mmm<0
                            fprintf('y{%d}(%d)=%1.5e<0\n',g,k,y{g}(k))
                            mmm = abs(mmm);
                        end
                        r = coef*(mmm)^(1/3)/2;
                        cubemass([coord(i,1), coord(i,2), coord(i,3)],r,'.b')
                        break
                    end
                end
            end
        end
    end
end

function cubemass(c,r,dotspec)
    delta = r/10;
    [X,Y,Z] = ndgrid(-r:delta:r,-r:delta:r,-r:delta:r) ; 
    plot3(c(1)+X(:),c(2)+Y(:),c(3)+Z(:),dotspec)
end