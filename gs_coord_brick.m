function coord = gs_coord_brick(nx,ny,nz,Lx,Ly,Lz)
    index = 1;
    coord = zeros(nx*ny*nz,3);
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                coord(index,:) = [i-1 j-1 k-1];
                index = index + 1;
            end
        end
    end
    coord(:,1) = coord(:,1)/max(coord(:,1))*Lx;
    coord(:,2) = coord(:,2)/max(coord(:,2))*Ly;
    coord(:,3) = coord(:,3)/max(coord(:,3))*Lz;
end