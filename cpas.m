function [x,A,b,exitflag] = cpas(c,A,b,lb,ub,nfix,data,verbose)
    
    report_iter = 0;
    if report_iter
        maxiter = 100000;
        fobj_iter = zeros(maxiter,1);
        % nrows_iter = zeros(maxiter,1);
        nrows2add_iter = cell(1,data.ncut);
        for k = 1:data.ncut
            nrows2add_iter{k} = zeros(maxiter,1);
        end
    end
    
    iterprint = 100;
    iterclear = 1000;
    cplex = Cplex('');
    cplex.Model.sense = 'minimize';
    cplex.addCols(c,[],lb,ub);
    cplex.addRows(-Inf*ones(length(b),1), A, b);
    A = cplex.Model.A;
    b = cplex.Model.rhs;
    cplex.DisplayFunc = @(s) ''; % suppress display output
    cplex.solve();
    exitflag = cplex.Solution.status;
    if exitflag > 1 && exitflag~=5
        fprintf('cplex.Solution.status %d.\n',exitflag)
        fprintf('%s.\n',cplex.Solution.statusstring)
        fprintf('Method = %d.\n',cplex.Solution.method)
        keyboard
    end
    x = cplex.Solution.x; % solves initial relaxation (x = argmin{c'*x : A*x<=b, lb<=x<=ub})
    nvar = length(x);
    iter = 1;


    % elapsedtimeid = tic;
    if verbose
        fprintf(['************************************************\n' ...
                 'CPAS\n' ...
                 '************************************************\n'])
        fprintf('Number of variables: %d\n',nvar)
        fprintf('Number of fix linear constraints: %d\n',nfix)
        fprintf('Number of cut functions: %d\n',data.ncut)
        fprintf('************************************************\n')
        fprintf('%d fobj=%1.9e nrows=%d\n',iter,c'*x,size(A,1))
    end
    % if os == "win"
    %     [user,~] = memory;
    %     bigmem = user.MemUsedMATLAB/2^20;
    % end

    while 1

        % itertimeid = tic;

        if mod(iter,iterclear) == 0
            % In truss structure A is dense since the eigenvalues are dense
            % cplex.Model.A is stored in sparse structure.
            % cplex.delRows do not dispose memory in cplex.Model.A
            % whos("cplex.Model.A").bytes >> numel(full(cplex.Model.A))*8.
            % To free memory in matlab:
            % 1. create variable fullA = full(cplex.Model.A)
            % 2. clear cplex structure
            % 3. restart cplex with fullA
            fullA = full(cplex.Model.A);
            % fprintf('fullA is %d bytes, clex.Model.A is %d bytes\n',numel(fullA)*8,whos("A").bytes)
            % if whos("A").bytes > 2^30
            %     fprintf('clex.Model.A is greater tha 1 GB !!!!!!!\n')
            %     keyboard
            % end
            clear cplex
            cplex = Cplex('');
            cplex.Model.sense = 'minimize';
            cplex.addCols(c,[],lb,ub);
            cplex.addRows(-Inf*ones(length(b),1), fullA, b); % if A instead of fullA, memory leaks
            A = cplex.Model.A;
            b = cplex.Model.rhs;
            cplex.DisplayFunc = @(s) ''; % suppress display output

            % if os == "win"
            %     [user,~] = memory;
            %     nowmem = user.MemUsedMATLAB/2^20;
            %     %fprintf('nowmem=%1.4e\n',nowmem)
            %     if nowmem>bigmem
            %        fprintf('MemUsedMATLAB increase %3.3e percent in 1000 iterations !!!\n',(nowmem-bigmem)/bigmem*100)
            %        bigmem = nowmem;
            %     end
            % end
            if strlength(data.file_mat_str) > 0
                % fprintf('Saving c,A,b,lb,ub,nfix,data.\n')
                save(data.file_mat_str,"c","fullA","b","lb","ub","nfix","data")
            end
        end
        
        cuts = [];
        for k = 1:data.ncut
            cuts = [cuts; data.cut{k}(x,data)];
            if report_iter
                nrows2add_iter{k}(iter) = size(cuts,1);
            end
        end
        
        if report_iter
            fobj_iter(iter) = c'*x;
            % nrows_iter(iter) = length(b);
        end

        if isempty(cuts)
            exitflag = 1;
            % elapsedtime = toc(elapsedtimeid);
            break
        end

        % if any(cuts(:,1:nvar)*x <= cuts(:,nvar+1))
        %     keyboard
        % end
        
        AA = cuts(:,1:nvar);
        bb = cuts(:,nvar+1);
        sa = sqrt(diag(AA*AA'));
        if verbose && mod(iter,iterprint)==0
            dist = min(abs(AA*x - bb)./sa);
            nrows = length(b);
            nrows2add = length(bb);
            fprintf('%d fobj=%1.9e nrows=%d nrows2add=%d dist=%1.9e \n', ...
                     iter,   c'*x, nrows,   nrows2add,   dist)
            % if os == "win"
            %     [user,~] = memory;
            %     matlabmem = user.MemUsedMATLAB/(1024*1024);
            %     fprintf('%d fobj=%1.9e dist=%1.9e nrows=%d nrows2add=%d matlabmem=%1.3eMB\n', ...
            %              iter,   c'*x,      dist,      nrows,   nrows2add,   matlabmem)
            % else
            %     fprintf('%d fobj=%1.9e dist=%1.9e nrows=%d nrows2add=%d\n', ...
            %              iter,   c'*x,      dist,      nrows,   nrows2add)
            % end
        end

        cplex.addRows(-Inf*ones(size(bb)), AA, bb);
        A = cplex.Model.A;
        b = cplex.Model.rhs;
        cplex.Start.x = x;
        cplex.solve();
        exitflag = cplex.Solution.status;
        if exitflag > 1 && exitflag~=5
            fprintf('cplex.Solution.status %d.\n',exitflag)
            fprintf('%s.\n',cplex.Solution.statusstring)
            fprintf('Method = %d.\n',cplex.Solution.method)
            keyboard
        end
        x = cplex.Solution.x;
        if size(A,1)-nfix > 4*nvar % 1200
            s = b - A*x - abs(cplex.Solution.dual);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
            s(1:nfix) = -1;
            [~,idx] = sort(s,'descend');
            numdel = size(A,1)-nfix-2*nvar; % - 900;
            whichdel = idx(1:numdel);
            idx_numdel = idx(numdel);
            least_del_slack = s(idx_numdel);
            % nbytes0 = whos("A").bytes;
            cplex.delRows(whichdel);
            A = cplex.Model.A;
            b = cplex.Model.rhs;
            % nbytes1 = whos("A").bytes;
            % if nbytes1 >= nbytes0
            %     fprintf('delRows do not release rows from memory !!!\n')
            % end
            cplex.solve();
            exitflag = cplex.Solution.status;
            if exitflag > 1 && exitflag~=5
                fprintf('cplex.Solution.status %d.\n',exitflag)
                fprintf('%s.\n',cplex.Solution.statusstring)
                fprintf('Method = %d.\n',cplex.Solution.method)
                keyboard
            end
            x0 = x;
            x = cplex.Solution.x;
            if norm(x-x0) > 1e-7
                fprintf('Remove %d rows and solution change %1.2e!\n',numdel,norm(x-x0))
                fprintf('fobj1-fobj0=%1.1e\n',c'*(x-x0))
                fprintf('least_del_slack = %1.6e\n',least_del_slack)
                % pause(1)
            end
        end

        iter = iter + 1;
        % elapsedtime = toc(elapsedtimeid);
        % if elapsedtime > 24*3600%60*5
        %     fprintf('CPAS: time exceeded... aborting\n')
        %     exitflag = -999;
        %     break
        % end
    end
    if verbose
        fprintf('%d iterations.\n',iter)
        % fprintf('Elapsed time: %d s.\n',elapsedtime)
        % fprintf('Iteration mean time: %d s.\n',itertimesum/iter)
        fprintf('CPAS returns %d cuts, exitflag %d.\n',length(b),exitflag)
        fprintf('************************************************\n')
    end

    if report_iter
        fobj_iter = fobj_iter(1:iter);
        % nrows_iter = nrows_iter(1:iter);
        for k = 1:data.ncut
            nrows2add_iter{k} = nrows2add_iter{k}(1:iter);
        end
        % save("report_iter.mat", "fobj_iter","nrows_iter","nrows2add_iter")
        save("report_iter.mat", "fobj_iter","nrows2add_iter")
    end
end