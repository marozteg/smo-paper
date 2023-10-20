function draw_truss(coord,elem,volume,show,tol,line_spec)
    if nargin == 3
        show = [1 1];
        tol = 1e-3;
        line_spec = '-k';
    elseif nargin == 4
        tol = 1e-3;
        line_spec = '-k';
    elseif nargin == 5
        line_spec = '-k';
    end
    nel = size(elem,1);
    gapx = 0.01*max(max(coord));
    gapy = gapx;
    gapz = gapy;
    fs = 11;
    if show(1)>0
        for i = 1:size(coord,1)
            plot3(coord(i,1), coord(i,2), coord(i,3),'.k')
            if show(1) == 2
                show_node_number(coord,i,gapx,gapy,fs)
            end
        end
    end
    if show(2)>0
        no1 = elem(:,1);
        no2 = elem(:,2);
        area = sqrt(volume./gs_length(coord,elem));
        for i = 1:nel
            cx = [coord(no1(i),1) coord(no2(i),1)];
            cy = [coord(no1(i),2) coord(no2(i),2)];
            cz = [coord(no1(i),3) coord(no2(i),3)];
            if volume(i) > tol
                % draw_elem_2d(coord(no1(i),:),coord(no2(i),:),area(i),line_spec)
                plot3(cx,cy,cz,line_spec,'linewidth',area(i))
            end
            if show(2) == 2
                show_elem_number(coord,no1(i),no2(i),gapx,gapy,fs,i)
            end
        end
    end
    xmin = min(coord(:,1));ymin = min(coord(:,2));zmin = min(coord(:,3));
    xmax = max(coord(:,1));ymax = max(coord(:,2));zmax = max(coord(:,3));
    xlim([xmin-gapx xmax+gapx])
    ylim([ymin-gapy ymax+gapy])
    zlim([zmin-gapz zmax+gapz])
    axis off%n
    axis equal
end

function show_node_number(coord,node,gx,gy,fs)
    text(coord(node,1)+gx,coord(node,2)+gy,coord(node,3),num2str(node),'FontSize',fs,'Color','m')
end

function show_elem_number(coord,n1,n2,gx,gy,fs,i)
    v = [coord(n2,1)-coord(n1,1),coord(n2,2)-coord(n1,2),coord(n2,3)-coord(n1,3)];
    lamb = 0.2;
    text(coord(n1,1)+lamb*v(1)+gx,coord(n1,2)+lamb*v(2)+gy,coord(n1,3)+lamb*v(3),num2str(i),'FontSize',fs,'Color','b')
end