clc
clear
close all
% % LINUX *************************************************************
%     addpath('/opt/ibm/ILOG/CPLEX_Studio1210/cplex/matlab/x86-64_linux')
%     run("/home/marozteg/MATLAB-Drive/fminsdp/install.m")
%     run("/home/marozteg/MATLAB-Drive/PENLABv104/install.m")
% % *******************************************************************
% WINDOWS ***********************************************************
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\matlab\x64_win64')
    run("C:\Users\maroz\MATLAB Drive\fminsdp\install.m")
    run("C:\Users\maroz\MATLAB Drive\PENLABv104\install.m")
% *******************************************************************

ver = 1;
linprog_options = optimoptions('linprog','Display','off');

% <GROUND_STRUCTURE>
for nnnx = 3:3 % 5:-1:3%3:5
    nnny = 2; %nnnx;
    nnnz = nnnx;
    LLLLx = 1.0;
    LLLLy = LLLLx/(nnnx-1);
    LLLLz = LLLLx;
    coord = gs_coord_brick(nnnx,nnny,nnnz,LLLLx,LLLLy,LLLLz);
    if nnnx == 2
        elem = gs_elem(coord,2,0);
    elseif nnnx == 3
        elem = gs_elem(coord,3,0); % ATT on overlap !!!!!
    elseif nnnx == 4
        elem = gs_elem(coord,6,0);
    elseif nnnx == 5
        elem = gs_elem(coord,3,0);
    elseif nnnx == 6
        elem = gs_elem(coord,11,0);
    end
    elem_pair = [];% elem_sym(coord,elem,1,3.0);
    L = gs_length(coord,elem);
    dof = gs_dof(elem);
    % figure;hold
    % draw_truss(coord,elem,ones(size(elem,1),1),[1 1],0,'k');
    % view(10,10)
% </GROUND_STRUCTURE>

% <SUPPORT>
    nn = size(coord,1);
    support_coord = [];
    kind = [];
    for i = 1:nn
        if coord(i,1) == 0 
            support_coord = [support_coord; coord(i,:)];
            kind = [kind; "xyz"];
        end
        % support_coord = [support_coord; coord(i,:)];
        % kind = [kind; "y"];
    end
    adof = support(coord, support_coord, kind);
% <\SUPPORT>

% <LOAD>
    kk = 1;
    clear load_case
    for i = 1:nn
        if (coord(i,1) == LLLLx) %&& (coord(i,2) == 0) % && (coord(i,3) == 0)
            % load_case{kk} = [coord(i,:) [-1 0 0]];
            % kk = kk + 1;
            load_case{kk} = [coord(i,:) [0 0 -1]];
            kk = kk + 1;
        end
    end
    F = loadcases(coord, load_case, adof);
% <\LOAD>

% <MATERIAL>
    E = 1.0;
    rho = 1.0;
% <\MATERIAL>

% <STIFFNESS_MATRIX>
    K0 = zeros(size(F,1),size(F,1));
    Kdx = mat_Kdx_vol(E*ones(size(elem,1),1),elem,coord,dof,adof);
    Kfun = @(x) sum_mat(1,size(elem,1),K0,Kdx,x);
% <\STIFFNESS_MATRIX>

% <MASS_MATRIX>
    Mdx = mat_Mdx_vol(rho*ones(size(elem,1),1),elem,coord,dof,adof);
    kk = 1;
    clear adofmass
    clear m0
    % MASS IN NODES TO THE RIGHT BOUNDARY
    % for i = 1:nn
    %     if (coord(i,1) == LLLLx)
    %         adofmass{kk} = zeros(3*nn,1);
    %         adofmass{kk}((i-1)*3+1:(i-1)*3+3) = 1;
    %         m0(kk) = 10;%*kk;
    %         kk = kk + 1;
    %     end
    % end
    % % MASS IN NODES TO THE BOTTON AND UPPER BOUNDARIES, EXCEPT LEFT AND RIGHT NODES
    % for i = 1:nn
    %     if (0 < coord(i,1) && coord(i,1) < LLLLx) && (coord(i,3) == 0 || coord(i,3) == LLLLz)
    %         adofmass{kk} = zeros(3*nn,1);
    %         adofmass{kk}((i-1)*3+1:(i-1)*3+3) = 1;
    %         m0(kk) = 10;%*kk;
    %         kk = kk + 1;
    %     end
    % end
    % MASS IN ALL NODES, EXCEPT SUPPORT NODES
    for i = 1:nn
        if coord(i,1) > 0 %&& coord(i,2) == 0 % && coord(i,3) == 0
            adofmass{kk} = zeros(3*nn,1);
            adofmass{kk}((i-1)*3+1:(i-1)*3+3) = 1;
            m0(kk) = 10;
            kk = kk + 1;
        end
    end
    nk = length(m0);
    M0 = cell(nk,1);
    Mfun = cell(nk,1);
    for k = 1:nk
        adofmass{k} = adofmass{k}(adof);
        M0{k} = diag(adofmass{k})*m0(k);
        Mfun{k} = @(x) sum_mat(1,size(elem,1),M0{k},Mdx,x);
    end
% <\MASS_MATRIX>

% <PROBLEM_SIZE>
    m = size(elem,1);
    [n,nlc] = size(F);
    nk = length(m0);
% <\PROBLEM_SIZE>

% <X0>
    area_0 = ones(m,1);
    x_0 = area_0.*L;
    vol_0 = sum(x_0);
    gamma_0 = zeros(nlc,1);
    for j = 1:nlc
        [u0,gamma_0(j)] = linprog(F(:,j),[],[],Kfun(x_0),F(:,j),[],[],linprog_options);
        %draw_truss(coord_desloc(coord,u0,adof),elem,x0,[0 1],0,'k','nofig');
    end
    [wcc_0,load_idx] = max(gamma_0);
    ss = strjoin(string(cell2mat(load_case(load_idx))));
    fprintf('worst case compliance is %1.4e and attained in load %s.\n',wcc_0,ss)
    lambda_0 = zeros(nk,1);
    for k = 1:nk
        lambda_0(k) = eigpair(Mfun{k}(x_0),Kfun(x_0),"min");    
    end
    [lambda_0,lamb_idx] = min(lambda_0);
    ss = strjoin(string(adof(adofmass{lamb_idx}==1)));
    fprintf('least eigenfreq is %1.4e at mass %d with dof %s.\n',lambda_0,lamb_idx,ss)
    figure;hold
    sz = min(L);
    draw_truss(coord,elem,x_0,[1 1],0,'k');
    draw_support(support_coord,kind,sz/5)
    draw_load(load_case,sz)
    m0_ = num2cell(m0');
    adofmass2 = cell(nk,1);    % draw_nonstructmass requires adofmass2, OMG
    for g = 1:nk
        adofmass2{g} = cell(1,1);
        adofmass2{g}{1} = adofmass{g};
    end
    draw_nonstructmass(coord,adofmass2,adof,m0_,sz)
    title(sprintf(['n=%d, m=%d, nl=%d, nk=%d\n' ...
                   'area=1, V=%1.4e, wcc=%1.4e, $\\lambda$=%1.4e\n'], ...
                   n,m,nlc,nk, ...
                   vol_0,max(gamma_0),lambda_0),'Interpreter','latex')
    view(10,10)
% <\X0>

% <DATA>
    data.vol = vol_0;
    data.gamma = wcc_0;
    data.f = F;
    data.K = Kfun;
    data.K0 = K0;
    data.Kdx = Kdx;
    data.M = Mfun;
    data.M0 = M0;
    data.Mdx = Mdx;
    data.file_mat_str = ""; % mfilename();
% <\DATA>


% <OPT_DATA>
    LAMBDA = 100;
    fprintf('max lambda\n')
    fprintf('x(%d),lambda(%d)\n',m,1)
    fprintf('s.t. [data.gamma -F(:,j)'';-F(:,j) K(x)]>=0, j=1,...,nlc(%d)\n',nlc)
    fprintf('     K(x)-lambda*(M(x)+M0^k)>=0, k=1,...,nk(%d)\n',nk)
    fprintf('     sum(x)<=data.vol(%1.3e) and\n',data.vol)
    fprintf('     lambda<=LAMBDA(%1.3e) and\n',LAMBDA)
    fprintf('     x>=0.\n')
    fprintf('n=%d, nnnx=%d, nnnz=%d, LLLLx=%1.3e, LLLLz=%1.3e\n',n,nnnx,nnnz,LLLLx,LLLLz)
    pause(3)
% <\OPT_DATA>

% <CPAS>
    c =  [zeros(m,1);-1];
    nep = size(elem_pair,1);
    Asym = zeros(nep,m);        % symmetry constraint
    for i = 1:nep
        Asym(i,elem_pair(i,1)) = 1;
        Asym(i,elem_pair(i,2)) = -1;
    end
    A =  [[ones(1,m) 0]
          [Asym zeros(nep,1)]];
    b =  [data.vol
          zeros(nep,1)];
    nfix = length(b); % constraints A(:,1:nfix)*x <= b(1:nfix) are "fixed": they will not be erased
    lb = [zeros(m,1);    -Inf];
    ub = [ones(m,1)*Inf; LAMBDA];  
    data.ncut = 2;
    data.cut{1} = @(x,data) compl_1(x,data);
    data.cut{2} = @(x,data) freq_1(x,data);
    data.TOL = 1e-6;
    data.options = cplexoptimset('cplex');

    tic
    [x_sol,A,b,exitflag_CPAS] = cpas(c,A,b,lb,ub,nfix,data,ver);
    % [x_sol,A,b,exitflag_CPAS] = cpas_load(mfilename(),1);
    time_CPAS = toc;

    if exitflag_CPAS > 0
        fprintf('CPAS time is %1.6e s\n',time_CPAS)
        x_CPAS = x_sol(1:m);
        vol_CPAS = sum(x_CPAS);
        area_CPAS = x_CPAS./L;
        gamma_CPAS = zeros(nlc,1);
        for j = 1:nlc
            [~,gamma_CPAS(j)] = linprog(F(:,j),[],[],Kfun(x_CPAS),F(:,j),[],[],linprog_options);
        end
        [wcc_CPAS,load_idx] = max(gamma_CPAS);
        ss = strjoin(string(cell2mat(load_case(load_idx))));
        fprintf('worst case compliance is %1.4e and attained in load %s.\n',wcc_CPAS,ss)
        lambda_CPAS = x_sol(m+1);
        lambda = zeros(nk,1);
        for k = 1:nk
            lambda(k) = eigpair(Mfun{k}(x_CPAS),Kfun(x_CPAS),"min");    
        end
        [lambda,lamb_idx] = min(lambda);
        ss = strjoin(string(adof(adofmass{lamb_idx}==1)));
        fprintf('least eigenfreq is %1.4e at mass %d with dof %s.\n',lambda,lamb_idx,ss)
        fprintf('lambda* is %1.4e\n',lambda_CPAS)
        figure;hold
        draw_truss(coord,elem,x_CPAS,[1 1],1e-3,'k');
        title(sprintf(['CPAS\n' ...
               'max $\\lambda$, n=%d, m=%d, nl=%d, nk=%d\n' ...
               '$sum(x^*)$=%1.4e, $wcc^*$=%1.4e, $\\lambda^*$=%1.4e\n' ...
               '$area_{min}^*$=%1.4e, $area_{max}^*$=%1.4e'], ...
               n,m,nlc,nk, ...
               vol_CPAS,wcc_CPAS,lambda_CPAS, ...
               min(area_CPAS),max(area_CPAS)),'Interpreter','latex')
        view(10,10)    
    else
        vol_CPAS = Inf;
        wcc_CPAS = Inf;
        lambda_CPAS = Inf;
        time_CPAS = Inf;
    end
% <\CPAS>


% <PENLAB>
    penm = [];
    penm.Nx = m+1;
    % objective function f
    penm.objfun  = @(x,Y,ud) deal([zeros(m,1);-1]'*x, ud);
    penm.objgrad = @(x,Y,ud) deal([zeros(m,1);-1],    ud);
    penm.objhess = @(x,Y,ud) deal(sparse(m+1,m+1),  ud);
    % constraint function g
    penm.NgLIN = 1; % volume constraint
    penm.confun  = @(x,Y,ud) deal([ones(m,1);0]'*x,ud);
    penm.congrad = @(x,Y,ud) deal([ones(m,1);0]',  ud);
    penm.ubg = data.vol;    % upper bound of g
    penm.lbg = -Inf;        % lower bound of g
    % matrix constraints A_{k}
    penm.NANLN = nk; % number of Bilinear Matrix Ineq constraints (BMI)
    penm.NALIN = nlc; % number of Linear Matrix Ineq constraints (LMI)
    penm.lbA = zeros(penm.NANLN+penm.NALIN,1); % A_i(x) must be positive semidef
    penm.mconfun  = @(x,Y,k,ud)     deal(mconfun( x,k,data),ud);
    penm.mcongrad = @(x,Y,k,i,ud)   deal(mcongrad(x,k,i,data),ud);
    penm.mconhess = @(x,Y,k,i,j,ud) deal(mconhess(x,k,i,j,data),ud);
    % x>=0
    penm.lbx = [zeros(m,1);-Inf]; % lower bound of x, -Inf gives no lang mult
    penm.ubx = [+Inf(m,1); LAMBDA];
    prob = penlab(penm);
    if ver == 1
        prob.opts = struct('outlev',2);
    else
        prob.opts = struct('outlev',0);
    end

    exitflag_PEN = prob.solve();

    if exitflag_PEN == 1
        time_PEN = prob.stats_time_total;
        fprintf('PENLAB time is %1.6e s\n',time_PEN)
        x_PEN = prob.x(1:m);
        vol_PEN = sum(x_PEN);
        area_PEN = x_PEN./L;
        gamma_PEN = zeros(nlc,1);
        for j = 1:nlc
            [~,gamma_PEN(j)] = linprog(F(:,j),[],[],Kfun(x_PEN),F(:,j),[],[],linprog_options);
        end
        [wcc_PEN,load_idx] = max(gamma_PEN);
        ss = strjoin(string(cell2mat(load_case(load_idx))));
        fprintf('worst case compliance is %1.4e and attained in load %s.\n',wcc_PEN,ss)
        lambda_PEN = prob.x(m+1);
        lambda = zeros(nk,1);
        for k = 1:nk
            lambda(k) = eigpair(Mfun{k}(x_CPAS),Kfun(x_CPAS),"min");    
        end
        [lambda,lamb_idx] = min(lambda);
        ss = strjoin(string(adof(adofmass{lamb_idx}==1)));
        fprintf('least eigenfreq is %1.4e at mass %d with dof %s.\n',lambda,lamb_idx,ss)
        fprintf('lambda* is %1.4e\n',lambda_PEN)
        figure;hold
        draw_truss(coord,elem,x_PEN,[1 1],1e-3,'k');
        title(sprintf(['PENLAB\n' ...
               'max $\\lambda$, n=%d, m=%d, nl=%d, nk=%d\n' ...
               '$sum(x^*)$=%1.4e, $wcc^*$=%1.4e, $\\lambda^*$=%1.4e\n' ...
               '$area_{min}^*$=%1.4e, $area_{max}^*$=%1.4e'], ...
               n,m,nlc,nk, ...
               vol_PEN,wcc_PEN,lambda_PEN, ...
               min(area_PEN),max(area_PEN)),'Interpreter','latex')
        view(10,10)
    else
        vol_PEN = Inf;
        wcc_PEN = Inf;
        lambda_PEN = Inf;
        time_PEN = Inf;
    end
% <\PENLAB>

% % <FMINSDP>
%     objfun = @(x) fminsdp_objfun(x);
%     A = [ones(1,m) 0];
%     b = data.vol;
%     lb = [zeros(m,1); 0];
%     ub = [+Inf(m,1); LAMBDA];
%     nonlcon = @(x) fminsdp_nonlcon(x,data);
%     if ver == 1
%         verb_str = 'iter-detailed';
%     else
%         verb_str = 'off';
%     end
%     svec_size = n*(n+1)/2;
%     Aind = zeros(1,nlc+1);
%     for i = 1:nlc+nk
%         Aind(i) = (i-1)*svec_size+1;
%     end
%     options = sdpoptionset(...
%         'Algorithm','interior-point',...
%         'GradObj','on', ...
%         'GradConstr','on', ...
%         'DerivativeCheck','off',...
%         'Display',verb_str, ...
%         'Aind',Aind);
% 
%     if nnnx < 7
%         tic
%         [x_sol,~,exitflag_FMS,~,~,~,~] = fminsdp(objfun,[x_0;lambda_0],A,b,[],[],lb,ub,nonlcon,options);
%         time_FMS = toc;
%     else
%         exitflag_FMS = -999;
%     end
% 
%     if exitflag_FMS == 1
%         fprintf('FMINSDP time is %1.6e s\n',time_FMS)
%         x_FMS = x_sol(1:m);
%         vol_FMS = sum(x_FMS);
%         area_FMS = x_FMS./L;
%         gamma_FMS = zeros(nlc,1);
%         for j = 1:nlc
%             [~,gamma_FMS(j)] = linprog(F(:,j),[],[],Kfun(x_FMS),F(:,j),[],[],linprog_options);
%         end
%         [wcc_FMS,load_idx] = max(gamma_FMS);
%         ss = strjoin(string(cell2mat(load_case(load_idx))));
%         fprintf('worst case compliance is %1.4e and attained in load %s.\n',wcc_FMS,ss)
%         lambda_FMS = x_sol(m+1);
%         lambda = zeros(nk,1);
%         for k = 1:nk
%             lambda(k) = eigpair(Mfun{k}(x_FMS),Kfun(x_FMS),"min");    
%         end
%         [lambda,lamb_idx] = min(lambda);
%         ss = strjoin(string(adof(adofmass{lamb_idx}==1)));
%         fprintf('least eigenfreq is %1.4e at mass %d with dof %s.\n',lambda,lamb_idx,ss)
%         fprintf('lambda* is %1.4e\n',lambda_FMS)
%         figure;hold
%         draw_truss(coord,elem,x_FMS,[1 1],1e-3,'k');
%         title(sprintf(['FMINSDP\n' ...
%                'max $\\lambda$, n=%d, m=%d, nl=%d, nk=%d\n' ...
%                '$sum(x^*)$=%1.4e, $wcc^*$=%1.4e, $\\lambda^*$=%1.4e\n' ...
%                '$area_{min}^*$=%1.4e, $area_{max}^*$=%1.4e'], ...
%                n,m,nlc,nk, ...
%                vol_FMS,wcc_FMS,lambda_FMS, ...
%                min(area_FMS),max(area_FMS)),'Interpreter','latex')
%         view(0,0)
%     else
%         vol_FMS = Inf;
%         wcc_FMS = Inf;
%         lambda_FMS = Inf;
%         time_FMS = Inf;
%     end
% % <\FMINSDP>

    % table(vol_0,vol_CPAS,vol_PEN,vol_FMS)
    % table(wcc_0,wcc_CPAS,wcc_PEN,wcc_FMS)
    % table(lambda_0,lambda_CPAS,lambda_PEN,lambda_FMS)
    % table(time_CPAS,time_PEN,time_FMS)


    table(vol_0,vol_CPAS,vol_PEN)%,vol_FMS)
    table(wcc_0,wcc_CPAS,wcc_PEN)%,wcc_FMS)
    table(lambda_0,lambda_CPAS,lambda_PEN)%,lambda_FMS)
    table(time_CPAS,time_PEN)%,time_FMS)

    save(strcat(mfilename(),"_nnnx_",num2str(nnnx)))

end


% <FMINSDP_HELP_FUNCTIONS>
function [fval, grad] = fminsdp_objfun(x)
    nvar = length(x);
    fval = - x(nvar);
    if nargout > 1
        grad = zeros(nvar,1);
        grad(nvar) = -1;
    end
end

function [cineq,ceq,cineqgrad,ceqgrad] = fminsdp_nonlcon(x,data)
    cineq = [];  % nonlinear inequality constraint
    nvar = length(x);
    F = data.f;
    [n,nlc] = size(F);
    nk = length(data.M0);
    svec_size = n*(n+1)/2;
    nmatcstr_svec = nlc*svec_size + nk*svec_size;
    ceq = zeros(nmatcstr_svec,1); % nonlinear equality constraint 
    for i = 1:nlc+nk
        if i <= nlc
            ceq((i-1)*svec_size+1 : i*svec_size) = svec(data.gamma*data.K(x) - F(:,i)*F(:,i)');
        elseif i > nlc
            k = i-nlc;
            ceq((i-1)*svec_size+1 : i*svec_size) = svec(data.K(x)-x(nvar)*data.M{k}(x));
        end
    end
    if nargout > 2
        cineqgrad = []; % nonlinear inequality constraint gradient
        ceqgrad = zeros(nvar,numel(ceq)); % nonlinear equality constraint gradient
        for i = 1:nvar-1
            Kdx = data.Kdx{i};
            Mdx = data.Mdx{i};
            ceqgrad(i, 1:nlc*svec_size) = kron(ones(1,nlc),svec(data.gamma*Kdx)');
            ceqgrad(i, nlc*svec_size+1:end) = kron(ones(1,nk),svec(Kdx - x(nvar)*Mdx)');
        end
        offset = nlc*svec_size;
        for k = 1:nk
            ceqgrad(nvar,offset+1+(k-1)*svec_size:offset+k*svec_size) = svec(-data.M{k}(x))';
        end
    end
end
% <\FMINSDP_HELP_FUNCTIONS>

% <PENLAB_HELP_FUNCTION>
function Ak = mconfun(x,k,data)
    nk = length(data.M0);
    if k <= nk
        Ak = sparse(data.K(x) - x(length(x))*data.M{k}(x)); % non LMI comes first
    elseif k > nk
        F = data.f;
        Ak = sparse([data.gamma -F(:,k-nk)'; -F(:,k-nk) data.K(x)]); % LMI
    end
end
function Aki = mcongrad(x,k,i,data)
    n = length(x); % [x;lambda]
    nk = length(data.M0);
    if k <= nk
        if i <= n - 1
            Aki = sparse(data.Kdx{i} - x(n)*data.Mdx{i}); % deriv of A_{k} wrt x(i)
        elseif i == n
            Aki = sparse(- data.M{k}(x)); % deriv of A_{k} wrt x(n)
        end
    elseif k > nk
        m = size(data.f,1);
        if i <= n - 1
            Aki = sparse([0 zeros(1,m); zeros(m,1) data.Kdx{i}]); % deriv of A_{k} wrt x(i)
        elseif i == n
            Aki = sparse(m+1,m+1); % deriv of A_{k} wrt x(n)
        end
    end
end
function Akij = mconhess(x,k,i,j,data)
    n = length(x);
    m = size(data.f,1);
    nk = length(data.M0);
    if k <= nk
        if i <= n - 1
            if j <= n - 1
                Akij = sparse(m,m);
            elseif j == n
                Akij = sparse(- data.Mdx{i});
            end
        elseif i == n
            if j <= n - 1
                Akij = sparse(- data.Mdx{j});
            elseif j == n
                Akij = sparse(m,m);
            end
        end
    elseif k > nk
        Akij = sparse(m+1,m+1);
    end
end
% <\PENLAB_HELP_FUNCTION>


% <CPAS_HELP_FUNCTION>
function allcuts = compl_1(x,data)
% maximize lambda
% x,lambda
% s.t.
%      data.gamma*K(x) -F(:,j)*F(:,j)'>=0, j=1,...,nlc   <<<------------
%      K(x)-lambda*M(x)>=0
%      sum(x)<=data.vol and
%      x>=0.
% uses one or more eigenvectors with negative eigenval
    F = data.f;
    nlc = size(F,2);
    allcuts = [];
    for j = 1:nlc
        [allEigVec,allEigVal] = eig(data.gamma*data.K(x)-F(:,j)*F(:,j)');
        mineig = min(diag(allEigVal));
        negIdx = find(diag(allEigVal) < min(-data.TOL,mineig/2));
        nneigenval = length(negIdx);
        if nneigenval == 0
            continue
        end
        W = allEigVec(:,negIdx);
        %W(abs(W)<1e-10) = 0;
        nvar = length(x);
        cuts = zeros(nneigenval, nvar+1);
        for ii = 1:nneigenval
            w = W(:,ii);
            for jj = 1:nvar-1
                cuts(ii,jj) = - data.gamma*w'*data.Kdx{jj}*w;
            end
            %cuts(i,nvar) = 0; % cuts(nvar)*vol
            cuts(ii,nvar+1) = data.gamma*w'*data.K0*w - (F(:,j)'*w)^2;
        end
        allcuts = [allcuts;cuts];
    end
end

function allcuts = compl_11(x,data)
% maximize lambda
% x,lambda
% s.t.
%      [data.gamma -F(:,j)';-F(:,j) K(x)]>=0, j=1,...,nlc   <<<------------
%      K(x)-lambda*M(x)>=0
%      sum(x)<=data.vol and
%      x>=0.
% uses one or more eigenvectors with negative eigenval
    allcuts = [];
    F = data.f;
    nlc = size(F,2);
    for j = 1:nlc
        [allEigVec,allEigVal] = eig([data.gamma -F(:,j)';-F(:,j) data.K(x)]);
        mineig = min(diag(allEigVal));
        negIdx = find(diag(allEigVal) < min(-data.TOL,mineig/2));
        nneigenval = length(negIdx);
        if nneigenval == 0
            continue
        end
        W = allEigVec(:,negIdx);
        %W(abs(W)<1e-10) = 0;
        nvar = length(x);
        cuts = zeros(nneigenval, nvar+1);
        for ii = 1:nneigenval
            h = W(1,ii);
            w = W(2:end,ii);
            for jj = 1:nvar-1
                cuts(ii,jj) = - w'*data.Kdx{jj}*w;
            end
            %cuts(i,nvar) = 0; % cuts(nvar)*lambda
            cuts(ii,nvar+1) = w'*data.K0*w + h*h*data.gamma - 2*h*F(:,j)'*w;
            cuts(i,nvar+1) = cuts(i,nvar+1) - data.TOL;
        end
        allcuts = [allcuts;cuts];
    end
end

function allcuts = freq_1(x,data) % uses all eigenvectors with negative eigenval
% maximize lambda
% x,lambda
% s.t.
%      data.gamma*K(x) -F(:,j)*F(:,j)'>=0, j=1,...,nl
%      K(x)-lambda*M^(k)(x)>=0             k=1,...,nk    <<<------------
%      sum(x)<=data.vol and
%      x>=0.
% uses one or more eigenvectors with negative eigenval
    allcuts = [];
    nvar = length(x);
    lambda = x(nvar);
    nk = size(data.M,1);
    for k = 1:nk
        [allEigVec,allEigVal] = eig(data.K(x) - lambda*data.M{k}(x));
        mineig = min(diag(allEigVal));
        negIdx = find(diag(allEigVal) < min(-data.TOL,mineig/2));
        nneigenval = length(negIdx);
        if nneigenval == 0
            continue
        end
        W = allEigVec(:,negIdx);
        % % w(abs(w)<1e-10) = 0;
        cuts = zeros(nneigenval, nvar+1);
        beta = zeros(nvar-1,1);
        for i = 1:nneigenval
            w0 = W(:,i);
            beta0 = w0'*data.M0{k}*w0;
            for j = 1:nvar-1
                beta(j) =  w0'*data.Mdx{j}*w0;
                cuts(i,j) = - (w0'*data.Kdx{j}*w0 - lambda*beta(j));
            end
            cuts(i,nvar) = max(beta)*data.vol + beta0; % cuts(nvar)*lambda
            cuts(i,nvar+1) = lambda*cuts(i,nvar) + w0'*data.K0*w0 - lambda*beta0;
            cuts(i,nvar+1) = cuts(i,nvar+1) - data.TOL;
        end
        allcuts = [allcuts;cuts];
    end
end
% <\CPAS_HELP_FUNCTION>