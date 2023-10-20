function coord = gs_coord_square(param,Lx)
    nx = 2*param + 1;
    coord = gs_coord_rectangle(nx,nx,Lx,Lx);
end