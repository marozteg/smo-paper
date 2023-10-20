function draw_support(support_coord,kind,height)
    color = 'g';
    lw = 2;
    for i = 1:size(support_coord,1)
        if contains(kind(i),"x")
            triline = xsupport(support_coord(i,:),height,height/2);
            for j = 1:3
                mat = triline{j};
                line(mat(1,:),mat(2,:),mat(3,:),'color',color,'linewidth',lw)
            end
        end
        if contains(kind(i),"y")
            triline = ysupport(support_coord(i,:),height,height/2);
            for j = 1:3
                mat = triline{j};
                line(mat(1,:),mat(2,:),mat(3,:),'color',color,'linewidth',lw)
            end
        end
        if contains(kind(i),"z")
            triline = zsupport(support_coord(i,:),height,height/2);
            for j = 1:3
                mat = triline{j};
                line(mat(1,:),mat(2,:),mat(3,:),'color',color,'linewidth',lw)
            end
        end
    end
end

function [triline] = xsupport(coord, height, base)
    prop = 0.5;
    triline = cell(3,1);
    triline{1} = [[coord(1) coord(1)-height];[coord(2) coord(2)];[coord(3) coord(3)+prop*base]];
    triline{2} = [[coord(1) coord(1)-height];[coord(2) coord(2)];[coord(3) coord(3)-prop*base]];
    triline{3} = [[coord(1)-height coord(1)-height];[coord(2) coord(2)];[coord(3)-prop*base coord(3)+prop*base]];
end
function [triline] = ysupport(coord, height, base)
    prop = 0.5;
    triline = cell(3,1);
    triline{1} = [[coord(1) coord(1)+prop*base];[coord(2) coord(2)+height];[coord(3) coord(3)]];
    triline{2} = [[coord(1) coord(1)-prop*base];[coord(2) coord(2)+height];[coord(3) coord(3)];];
    triline{3} = [[coord(1)-prop*base coord(1)+prop*base];[coord(2)+height coord(2)+height];[coord(3) coord(3)];];
end
function [triline] = zsupport(coord, height, base)
    prop = 0.5;
    triline = cell(3,1);
    triline{1} = [[coord(1) coord(1)+prop*base];[coord(2) coord(2)];[coord(3) coord(3)-height]];
    triline{2} = [[coord(1) coord(1)-prop*base];[coord(2) coord(2)];[coord(3) coord(3)-height]];
    triline{3} = [[coord(1)-prop*base coord(1)+prop*base];[coord(2) coord(2)];[coord(3)-height coord(3)-height]];
end