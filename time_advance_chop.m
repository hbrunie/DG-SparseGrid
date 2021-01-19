function f = time_advance_chop(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax)

if strcmp(opts.timestep_method,'BE')
    % Backward Euler (BE) first order
    f = backward_euler(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'FE')
    % Forward Euler (FE) first order
    f = forward_euler(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'matrix_exponential')
    % Matrix Exponential (ME) all order
    f = matrix_exponential(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'time_independent')
    % time independent d/dt==0
    f = time_independent(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'CN')
    % Crank Nicolson (CN) second order
    f = crank_nicolson(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif sum(strcmp(opts.timestep_method,{'ode15i','ode15s','ode45','ode23s'}))>0
    % One of the atlab ODE integrators.
    f = ODEm(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
else
    % RK3 TVD
    f = RungeKutta3(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
end

end

%% Use Matlab ODE solvers

function fval = ODEm(pde,opts,A_data,f0,t0,dt,deg,hash_table,Vmax,Emax)

clear strCR;

applyLHS = ~isempty(pde.termsLHS);
if applyLHS
    error('ERROR: Matlab ODE integrators not yet implemented for LHS=true');
end

%     function res = fun_for_jacobianest(x)
%     res = explicit_ode(t0,x);
%     disp('calling ode');
%     end

    function dfdt = explicit_ode(t,f)
        bc = boundary_condition_vector(pde,opts,hash_table,t);
        source = source_vector(pde,opts,hash_table,t);
        f1 = apply_A(pde,opts,A_data,f,deg,Vmax,Emax);
        dfdt = f1 + source + bc;
    end

    function res = implicit_ode(t,f,dfdt)
        bc = boundary_condition_vector(pde,opts,hash_table,t);
        source = source_vector(pde,opts,hash_table,t);
        f1 = apply_A(pde,opts,A_data,f,deg,Vmax,Emax);
        res = dfdt - (f1 + source + bc);
    end

if strcmp(opts.timestep_method,'ode45')
    
    stats = 'off';
    if(~opts.quiet)
        disp('Using ode45');
        stats = 'on';
    end
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'Stats',stats);
    [tout,fout] = ode45(@explicit_ode,[t0 t0+dt],f0,options);
    
elseif strcmp(opts.timestep_method,'ode15s')
    
    %     % estimate Jacobian
    %     numjacopts.diffvar = 2;
    %     numjacopts.vectvars = [];
    %     numjacopts.thresh = 1e-10;
    %     numjacopts.fac = [];
    
    %     disp('running odenumjac')
    %     rhs = feval(@explicit_ode,t0,f0);
    %     J = odenumjac(@explicit_ode,{t0 f0},rhs,numjacopts);
    %     S = sparse(J~=0.0);
    %     disp('done')
    
    %     disp('running jacobianest')
    %     [J2,err] = jacobianest(@fun_for_jacobianest,f0);
    %     disp('done');
    
    % call ode15s
    stats = 'off';
    output_func = '';   
    if(~opts.quiet)
        disp('Using ode15s');
        stats = 'on';
        output_func = @odetpbar;       
    end
    options = odeset('RelTol',1e-6,'AbsTol',1e-8,...
        'Stats',stats,'OutputFcn',output_func,'Refine',20);%,'Jacobian', J2);%'JPattern',S);
    [tout,fout] = ode15s(@explicit_ode,[t0 t0+dt],f0,options);

elseif strcmp(opts.timestep_method,'ode23s')
    
    % call ode23s
    stats = 'off';
    output_func = '';   
    if(~opts.quiet)
        disp('Using ode23s');
        stats = 'on';
        output_func = @odetpbar;       
    end
    options = odeset('Stats',stats,'OutputFcn',output_func);
    [tout,fout] = ode23s(@explicit_ode,[t0 t0+dt],f0,options);
    
elseif strcmp(opts.timestep_method,'ode15i')
    
    dfdt0 = f0.*0;
    [f0,dfdt0,resnrm] = decic(@implicit_ode,t0,f0,f0.*0+1,dfdt0,[]);
    if(~opts.quiet);disp('Using ode15i');end
    [tout,fout] = ode15i(@implicit_ode,[t0 t0+dt],f0,dfdt0);
    
end

fval = reshape(fout(end,:),size(f0));

end

%% 3-rd Order Kutta Method (explicit time advance)
function fval = RungeKutta3(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax)
options.format = 'h'; options.round = 1; options.subnormal = 1;
chop([],options)
%Convert pde to single
%disp(fieldnames(pde));
for d=1:numel(pde.dimensions)
    pde_s.dimensions{d}.name = pde.dimensions{d}.name;
    pde_s.dimensions{d}.min = pde.dimensions{d}.min;
    pde_s.dimensions{d}.max = pde.dimensions{d}.max;
    pde_s.dimensions{d}.lev = pde.dimensions{d}.lev;
    %function handle
    %pde_s.dimensions{d}.init_cond_fn =pde.dimensions{d}.init_cond_fn;
    %pde_s.dimensions{d}.jacobian = chop(pde.dimensions{d}.jacobian);
    pde_s.dimensions{d}.init_cond_fn =pde.dimensions{d}.init_cond_fn;
    pde_s.dimensions{d}.jacobian =pde.dimensions{d}.jacobian;
end

%COnversion from MD_TERM to chop not possible
%for d=1:numel(pde.terms)
%    pde_s.terms{d}      = chop(pde.terms{d});
%end
pde_s.terms      = pde.terms;

%Conversion to chop from struct is not possible
%for d=1:numel(pde.params)
%    pde_s.params     = chop(pde.params);
%end
pde_s.params     =pde.params;

for d=1:numel(pde.sources)
    pde_s.sources{d}    = chop(pde.sources{d});
end
pde_s.sources    =pde.sources;
pde_s.transform_blocks    =pde.transform_blocks;

pde_s.termsLHS   = chop(pde.termsLHS);
for d=1:numel(pde.boundary_conditions)
    pde_s.boundary_conditions{d} = chop(pde.boundary_conditions{d});
end
%Conversion to chop from cell is not possible
%for d=1:numel(pde.initial_conditions)
%    pde_s.initial_conditions{d} = chop(pde.initial_conditions{d});
%end
%pde_s.solutions = chop(pde.solutions);
%pde_s.connectivity = chop(pde.connectivity);
pde_s.initial_conditions =pde.initial_conditions;
pde_s.solutions =pde.solutions;
pde_s.connectivity =pde.connectivity;
num_terms     = numel(pde.terms);
num_terms_LHS = numel(pde.termsLHS);
num_dims      = numel(pde.dimensions);
for t=1:num_terms
    for d=1:num_dims
        pde_s.terms{t}.terms_1D{d}.mat = chop(pde.terms{t}.terms_1D{d}.mat);
    end
end


%Convert A_data to chop
dims = pde.dimensions;
nDims = numel(dims);
A_data_s.element_global_row_index = chop(A_data.element_global_row_index);
for d=1:nDims
    A_data_s.element_local_index_D{d} = chop(A_data.element_local_index_D{d});
end
%convert f, t, dt, deg, hash_table, Vmax,Emax to chop
f_s = chop(f);
t_s = chop(t);
dt_s = chop(dt);
deg_s = chop(deg);
%struct not possible
%hash_table_s = chop(hash_table);
hash_table_s =hash_table;
Vmax_s = chop(Vmax);
Emax_s = chop(Emax);

%disp(sprintf('RK3 class of A_data is %s', class(A_data.element_global_row_index) ));
%for d=1:nDims
%    disp(sprintf('RK3 class of A_data is %s', class(A_data.element_local_index_D{d}) ));
%end
%disp(sprintf('RK3 class of A_data_s is %s', class(A_data_s.element_global_row_index) ));
%for d=1:nDims
%    disp(sprintf('RK3 class of A_data_s is %s', class(A_data_s.element_local_index_D{d}) ));
%end

%%
% Sources
c2 = chop(1/2); c3 = chop(1);
source1 = source_vector(pde_s,opts,hash_table_s,t_s);
source2 = source_vector(pde_s,opts,hash_table_s,t_s+c2*dt_s);
source3 = source_vector(pde_s,opts,hash_table_s,t_s+c3*dt_s);

%%
% Inhomogeneous dirichlet boundary conditions
bc1 = boundary_condition_vector(pde_s,opts,hash_table_s,t_s);
bc2 = boundary_condition_vector(pde_s,opts,hash_table_s,t_s+c2*dt_s);
bc3 = boundary_condition_vector(pde_s,opts,hash_table_s,t_s+c3*dt_s);

% %%
% Apply any non-identity LHS mass matrix coefficient

applyLHS = ~isempty(pde.termsLHS);

a21 = chop(1/2); a31 = chop(-1); a32 = chop(2);
b1 = chop(1/6); b2 =  chop(2/3); b3 = chop(1/6);
%disp_x = [a21,a31,a32,b1,b2,b3]
%disp(disp_x)

if applyLHS
    [k1,A1,ALHS] = apply_A_chop(pde_s,opts,A_data_s,f_s,deg_s,Vmax_s,Emax_s);
    rhs1 = source1 + bc1;
    %     invMatLHS = inv(ALHS); % NOTE : assume time independent for now for speed.
    %     k1 = invMatLHS * (k1 + rhs1);
    k1 = ALHS \ (k1 + rhs1);
    y1 = f_s + dt_s*a21*k1;
    
    [k2] = apply_A_chop(pde_s,opts,A_data_s,y1,deg_s,Vmax_s,Emax_s);
    rhs2 = source2 + bc2;
    %     k2 = invMatLHS * (k2 + rhs2);
    k2 = ALHS \ (k2 + rhs2);
    y2 = f_s + dt_s*(a31*k1 + a32*k2);
    
    k3 = apply_A_chop(pde_s,opts,A_data_s,y2,deg_s,Vmax_s,Emax_s);
    rhs3 = source3 + bc3;
    %     k3 = invMatLHS * (k3 + rhs3);
    k3 = ALHS \ (k3 + rhs3);
else
    k1 = apply_A_chop(pde_s,opts,A_data_s,f_s,deg_s,Vmax_s,Emax_s)  + source1 + bc1;
    y1 = f_s + dt_s*a21*k1;
    k2 = apply_A_chop(pde_s,opts,A_data_s,y1,deg_s,Vmax_s,Emax_s) + source2 + bc2;
    y2 = f_s + dt_s*(a31*k1 + a32*k2);
    k3 = apply_A_chop(pde_s,opts,A_data_s,y2,deg_s,Vmax_s,Emax_s) + source3 + bc3;
end

fval = double(f_s + dt_s*(b1*k1 + b2*k2 + b3*k3));
%disp(sprintf('RK3 class of all data is %s %s %s %s %s %s %s %s', class(f),class(dt),class(k1),class(b2),class(k2),class(k3),class(b3),class(fval) ));
%
%disp(sprintf('RK3 class of A_data is %s', class(A_data.element_global_row_index) ));
%for d=1:nDims
%    disp(sprintf('RK3 class of A_data is %s', class(A_data.element_local_index_D{d}) ));
%end
%disp(sprintf('RK3 class of A_data_s is %s', class(A_data_s.element_global_row_index) ));
%for d=1:nDims
%    disp(sprintf('RK3 class of A_data_s is %s', class(A_data_s.element_local_index_D{d}) ));
%end

end

%% Time independent solve d/dt==0
function f1 = time_independent(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

s0 = source_vector(pde,opts,hash_table,t+dt);
bc0 = boundary_condition_vector(pde,opts,hash_table,t+dt);

[~,A] = apply_A(pde,opts,A_data,f0,deg);

f1 = -A \ (s0+bc0);

end

%% Matrix Exponential
function f1 = matrix_exponential(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

% Note that this is not implemented for source terms, so is only here for
% testing purposed. Adding sources terms requires adding another term,
% i.e., 
% f1 = expm(t*A) * int_0^t exp(-u*A) s(t) du  + expm(t*A)*f0
% as well as something for the boundary terms - that integral looks like a
% pain.

s0 = source_vector(pde,opts,hash_table,t+dt);
bc0 = boundary_condition_vector(pde,opts,hash_table,t+dt);

applyLHS = ~isempty(pde.termsLHS);

if applyLHS
    error('apply LHS not implemented for FE');
else
    [~,A] = apply_A(pde,opts,A_data,f0,deg);
    f1 = expm(A*dt)*f0; 
end

end

%% Forward Euler
function f1 = forward_euler(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

s0 = source_vector(pde,opts,hash_table,t+dt);
bc0 = boundary_condition_vector(pde,opts,hash_table,t+dt);

applyLHS = ~isempty(pde.termsLHS);

if applyLHS
    error('apply LHS not implemented for FE');
else
    f1 = f0 + dt * (apply_A(pde,opts,A_data,f0,deg) + s0 + bc0);
end

end

%% Backward Euler (first order implicit time advance)
function f1 = backward_euler(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

applyLHS = ~isempty(pde.termsLHS);

s0 = source_vector(pde,opts,hash_table,t+dt);
bc0 = boundary_condition_vector(pde,opts,hash_table,t+dt);

if opts.time_independent_A || opts.time_independent_build_A
    persistent A;
end

if opts.time_independent_A % relies on inv() so is no good for poorly conditioned problems
    persistent AA_inv;
    persistent ALHS_inv;
    
    if isempty(AA_inv)
        if applyLHS
            [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
            I = speye(numel(diag(A)));
            ALHS_inv = inv(ALHS);
            AA = I - dt*(ALHS_inv * A);
            b = f0 + dt*(ALHS_inv * (s0 + bc0));
        else
            [~,A] = apply_A(pde,opts,A_data,f0,deg);
            I = speye(numel(diag(A)));
            AA = I - dt*A;
            b = f0 + dt*(s0 + bc0);
        end
        
        if numel(AA(:,1)) <= 4096
            condAA = condest(AA);
            if ~opts.quiet; disp(['    condest(AA) : ', num2str(condAA,'%.1e')]); end
            
            if condAA > 1e6
                disp(['WARNING: Using time_independent_A=true for poorly conditioned system not recommended']);
                disp(['WARNING: cond(A) = ', num2str(condAA,'%.1e')]);
            end
        end
        
        AA_inv = inv(AA);
        f1 = AA_inv * b;
    else
        if applyLHS
            b = f0 + dt*(ALHS_inv * (s0 + bc0));
        else
            b = f0 + dt*(s0 + bc0);
        end
        f1 = AA_inv * b;
    end
   
else % use the backslash operator instead
    if applyLHS
        if opts.time_independent_build_A 
            if isempty(A);[~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);end
        else
            [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
        end
        I = speye(numel(diag(A)));
        AA = I - dt*(ALHS \ A);
        b = f0 + dt*(ALHS \ (s0 + bc0));
    else
        if opts.time_independent_build_A
            if isempty(A);[~,A] = apply_A(pde,opts,A_data,f0,deg);end
        else
            [~,A] = apply_A(pde,opts,A_data,f0,deg);
        end
        I = speye(numel(diag(A)));
        AA = I - dt*A;
        b = f0 + dt*(s0 + bc0);
    end
    
    % Direct solve
    rescale = false;
    if rescale
        %     [AA_rescaled,diag_scaling] = rescale2(AA);
        %     AA_thresholded = sparsify(AA_rescaled,1e-5);
        [P,R,C] = equilibrate(AA);
        AA_rescaled = R*P*AA*C;
        b_rescaled = R*P*b;
        f1 = AA_rescaled \ b_rescaled;
        f1 = C*f1;
    else
        f1 = AA \ b;
    end
    
%     % Pre-kron solve
%     num_terms = numel(pde.terms);
%     num_dims = numel(pde.dimensions);
%     for d=1:num_dims
%         dim_mat_list{d} = zeros(size(pde.terms{1}.terms_1D{1}.mat));
%         for t=1:num_terms           
%             dim_mat_list{d} = dim_mat_list{d} + pde.terms{t}.terms_1D{d}.mat;
%         end
%     end
%     for d=1:num_dims
%         A = dim_mat_list{d};
%         I = speye(numel(diag(A)));
%         AA_d =  I - dt*A;
%         dim_mat_inv_list{d} = inv(AA_d);
%     end
%     use_kronmultd = true;
%     if use_kronmultd
%         f1a = kron_multd(num_dims,dim_mat_inv_list,b);
%     else
%         f1a = kron_multd_full(num_dims,dim_mat_inv_list,b);
%     end
    
%     % Iterative solve
%     restart = [];
%     tol=1e-6;
%     maxit=1000;
%     tic;
%     [f10,flag,relres,iter,resvec] = gmres(AA,b,restart,tol,maxit);
%     t_gmres = toc;
%     figure(67)
%     semilogy(resvec);
%   
%     % Direct solve - LU approach to reducing repeated work   
%     [L,U,P] = lu(AA);
%     n=size(AA,1);
%     ip = P*reshape(1:n,n,1); % generate permutation vector
%     err = norm(AA(ip,:)-L*U,1); % err should be small
%     tol = 1e-9;
%     isok = (err <= tol * norm(AA,1) ); % just a check
%     disp(['isok for LU: ', num2str(isok)]);
%     % to solve   a linear system    A * x = b, we have   P * A * x = P*b
%     % then from LU factorization, we have (L * U) * x = (P*b),  so    x  = U \ (L \ ( P*b))
%     f1_LU =   U \ (L \  (P*b));
%     
%     % Direct solve - QR reduced rank
%     n=size(AA,1);
%     [Q,R,P] = qr(AA); % A*P = Q*R
%     % where R is upper triangular,   Q is orthogonal, Q’*Q is identity, P is column permutation
%     err = norm( AA*P - Q*R,1);
%     tol=1e-9;
%     isok = (err <= tol * norm(AA,1));  % just a check
%     disp(['isok for QR: ', num2str(isok)]);
%     % to solve   A * x = b,   we have A * P * (P’*x) = b, (Q*R) * (P’*x) = b
%     % y =   R\(Q’*b),   P*y = x
%     tol=1e-9;
%     is_illcond = abs(R(n,n)) <= tol * abs(R(1,1));
%     if(is_illcond)
%         disp('is_illcond == true');
%         bhat = Q'*b;
%         y = zeros(n,1);
%         yhat =  R(1:(n-1), 1:(n-1)) \ bhat(1:(n-1));  % solve with smaller system
%         y(1:(n-1)) = yhat(1:(n-1));
%         y(n) = 0;  % force last component to be zero
%         f1_QR = P * y;
%     else
%         disp('is_illcond == false');
%         f1_QR = P * (R \ (Q'*b));
%     end
%       
%     disp(['f1-f10:  ',num2str(norm(f1-f10)/norm(f1))]);
%     disp(['f1-f1_LU:  ',num2str(norm(f1-f1_LU)/norm(f1))]);
%     disp(['f1-f1_QR:  ',num2str(norm(f1-f1_QR)/norm(f1))]);
%     disp(['direct runtime: ', num2str(t_direct)]);
%     disp(['gmres runtime: ', num2str(t_gmres)]);
%     
%     f1 = f1_QR;
    
end
end


%% Crank Nicolson (second order implicit time advance)
function f1 = crank_nicolson(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

applyLHS = ~isempty(pde.termsLHS);

s0 = source_vector(pde,opts,hash_table,t);
s1 = source_vector(pde,opts,hash_table,t+dt);

bc0 = boundary_condition_vector(pde,opts,hash_table,t);
bc1 = boundary_condition_vector(pde,opts,hash_table,t+dt);

if opts.time_independent_A % uses inv() so no good for poorly conditioned systems
    persistent A;
    persistent AA_inv;
    persistent ALHS_inv;
    
    if isempty(AA_inv)
        if applyLHS
            [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
            I = speye(numel(diag(A)));
            ALHS_inv = inv(ALHS);
            AA = 2*I - dt*(ALHS_inv * A);
            b = 2*f0 + dt*(ALHS_inv * A)*f0 + dt*(ALHS_inv * (s0+s1+bc0+bc1));
        else
            [~,A] = apply_A(pde,opts,A_data,f0,deg);
            I = speye(numel(diag(A)));
            AA = 2*I - dt*A;
            b = 2*f0 + dt*A*f0 + dt*(s0+s1) + dt*(bc0+bc1);
        end
        
        rcondAA = rcond(AA);
        if ~opts.quiet; disp(['    rcond(AA) : ', num2str(rcondAA)]); end
        
        if 1/rcondAA > 1e6
            disp(['WARNING: Using time_independent_A=true for poorly conditioned system not recommended']);
            disp(['WARNING: cond(A) = ', num2str(rcondAA)]);
        end
        
        AA_inv = inv(AA);
        f1 = AA_inv * b;
    else        
        if applyLHS
            b = 2*f0 + dt*(ALHS_inv * A)*f0 + dt*(ALHS_inv * (s0+s1+bc0+bc1));
        else
            b = 2*f0 + dt*A*f0 + dt*(s0+s1) + dt*(bc0+bc1);
        end    
        f1 = AA_inv * b;
    end
    
else % use the backslash operator for time_indepent_A = false
    if applyLHS
        [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
        I = speye(numel(diag(A)));
        AA = 2*I - dt*(ALHS \ A);
        b = 2*f0 + dt*(ALHS \ A)*f0 + dt*(ALHS \ (s0+s1+bc0+bc1));
    else
        [~,A] = apply_A(pde,opts,A_data,f0,deg);
        I = speye(numel(diag(A)));
        AA = 2*I - dt*A;
        b = 2*f0 + dt*A*f0 + dt*(s0+s1) + dt*(bc0+bc1);
    end
    f1 = AA \ b;
end

end
