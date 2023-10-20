function elem = gs_elem(coord,sz,overlap)

    if nargin == 2
        overlap = 1;
    end

    nn = size(coord,1);

    % all [a b] with a and b in 1:nn 
    [n1,n2] = meshgrid(1:nn,1:nn);
    elem = [n1(:) n2(:)];
    
    % Remove [a a]
    elem(elem(:,1)==elem(:,2),:) = [];
    
    % Remove [a b] if a>b
    elem = sort(elem,2);
    elem = unique(elem,'rows');
  
    if not(overlap)
        % Remove elements longer and in the same line than smaller ones
        nel = size(elem,1);
        elem2rem = false(nel,1);
        for i = 1:nel
            element = elem(i,:);
            n1 = element(1);
            n2 = element(2);
            v = [coord(n2,1)-coord(n1,1)
                 coord(n2,2)-coord(n1,2)
                 coord(n2,3)-coord(n1,3)];
            for n3 = 1:nn
                if n3 ~= n1 && n3 ~= n2
                    w = [coord(n3,1)-coord(n1,1)
                         coord(n3,2)-coord(n1,2)
                         coord(n3,3)-coord(n1,3)];
                    if norm(cross(v,w)) < 1e-6 % v and w are parallel
                        [~,j] = max(abs(v));
                        if 0 < w(j)/v(j) && w(j)/v(j) < 1
                            elem2rem(i) = true;
                        end
                    end
                end
            end
        end
        elem(elem2rem,:) = [];
    end

    % Remove bars longer than distmax
    n1 = elem(:,1);
    n2 = elem(:,2);
    c1x = coord(n1,1);
    c1y = coord(n1,2);
    c1z = coord(n1,3);
    c2x = coord(n2,1);
    c2y = coord(n2,2);
    c2z = coord(n2,3);
    dx = c2x-c1x;
    dy = c2y-c1y;
    dz = c2z-c1z;
    lengths = sqrt(dx.^2 + dy.^2 + dz.^2);

    sort_length = sort(lengths);
    % unique_sort_length = uniquetol(sort_length,1e-6);
    unique_sort_length = unique(sort_length);
    % % problem: due to roundoff there are many "unique" lengths
    % % todo: must cluster in groups of similar lengths
    % npoints = 100000;
    % hh = max(sort_length)/npoints;
    % group = 1;
    % xx = 0;
    % unique_sort_length = [];
    % for i = 1:npoints
    %     found = 0;
    %     idx_cluster = [];
    %     for j = 1:length(sort_length)
    %         if abs(sort_length(j)-xx)<hh
    %             idx_cluster = [idx_cluster j];
    %             found = 1;
    %         end
    %     end
    %     if ~found
    %         idx_cluster = unique(idx_cluster);
    %         unique_sort_length(group) = mean(sort_length(idx_cluster)); 
    %         group = group + 1;
    %     end
    %     xx = xx + hh;
    % end
    
    if sz > length(unique_sort_length)
        sz = length(unique_sort_length);
        fprintf('gs_elem: taking argument equal to %d.\n',sz)
    end

    distmax = unique_sort_length(sz);
    elem(lengths>distmax,:) = [];
end