function draw_load(load_cases,sz)
    color = 'r';
    prop = 0.35;
    lw = 4;
    nl = size(load_cases,2);
    for lc = 1:nl
        loadcase = cell2mat(load_cases(lc));
        load_coords = loadcase(:,1:3);
        vectors = loadcase(:,4:6);
        vectors = vectors./vecnorm(vectors')';
        vectors = sz*prop*vectors;
        for i = 1:size(loadcase,1)
            vect = vectors(i,:);
            foot = load_coords(i,:);
            head = foot + vect;
            orth_vect = [-vect(3) vect(2) vect(1)];
            line( ...
                [foot(1) head(1)], ...
                [foot(2) head(2)], ...
                [foot(3) head(3)], ...
                'color',color, ...
                'linewidth', lw ...
            );
            line( ...
                [head(1) head(1)-prop*(vect(1)+orth_vect(1))], ...
                [head(2) head(2)-prop*(vect(2)+orth_vect(2))], ...
                [head(3) head(3)-prop*(vect(3)+orth_vect(3))], ...
                'color',color, ...
                'linewidth', lw/2 ...
            );
            line( ...
                [head(1) head(1)-prop*(vect(1)-orth_vect(1))], ...
                [head(2) head(2)-prop*(vect(2)-orth_vect(2))], ...
                [head(3) head(3)-prop*(vect(3)-orth_vect(3))], ...
                'color',color, ...
                'linewidth', lw/2 ...
            );
        end
    end
end