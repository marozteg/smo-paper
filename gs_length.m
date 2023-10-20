function L = gs_length(coord,elem)
    L = vecnorm((coord(elem(:,2),:)-coord(elem(:,1),:))')';
end