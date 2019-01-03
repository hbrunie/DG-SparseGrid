function [err,fval,fval_realspace] = fk6d(pde,Lev,Deg,TEND,quiet,compression,implicit,gridType)

%% MATLAB (reference) version of the DG-SG solver
% [err,fval,fval_realspace] = fk6d(pde,Lev,Deg,TEND,quiet,compression,implicit,gridType)
%
% The default execution solves the Vlasov-Poisson system of equations
%
% $$f_t + v\frac{\partial f}{\partial x} + E\left(x,t\right)\frac{\partial
% f}{\partial v} f=0$$
%
% $$-\frac{\partial^2\phi}{\partial x^2} = \rho - 1$$
%
% $$E\left(x,t\right)= \frac{\partial\phi}{\partial x}$$
%
% where $\rho\left(x,t\right)=\int_v f(x,v,t) dv$.

format short e
folder = fileparts(which(mfilename));
addpath(genpath(folder));

%% Step 1. Set input parameters
% pde :  A structure containing the initial condition and domain
% information. See PDE/vlasov4.m for and example. Note that this does not
% actually describe the PDE (that's done by the coefficient matricies), so
% we should probably rename this.
%
% Lev: Maximum level of the mesh in all dimensions.
%
% Deg: Degree of basis functions.
%
% TEND : End time of the simulation, i.e., run from t=0 to t=TEND in steps
% of dt.
%
% quiet : Print debugging statements or not.
%
% Compression : Choice of approach to constructing the A system matrix.

if ~exist('pde','var') || isempty(pde)
    % Equation setup
    pde = Vlasov8;
end

nTerms = numel(pde.terms);
nDims = numel(pde.dimensions);

if ~exist('TEND','var') || isempty(TEND)
    % End time
    TEND = 0.001;
end
if exist('Lev','var')
    % Number of levels
    for d=1:nDims
        pde.dimensions{d}.lev = Lev;
    end
end

if exist('Deg','var')
    % Polynomial degree
    % Deg = 2 Means Linear Element
    for d=1:nDims
        pde.dimensions{d}.deg = Deg;
    end
end

if ~exist('quiet','var') || isempty(quiet)
    % Enable / disable print statements
    quiet = 0;
end
if ~exist('compression','var') || isempty(compression)
    % Use or not the compression reference version
    compression = 3;
end
if ~exist('gridType','var') || isempty(gridType)
    gridType = 'SG';%'FG'
else
    if strcmp(gridType,'SG') || strcmp(gridType,'FG')
    else
        error("gridType must be set to 'SG' or 'FG'");
    end
end
if ~exist('implicit','var') || isempty(implicit)
    pde.implicit = 0;
else
    pde.implicit = implicit;
end

if pde.implicit
    compression = 1;
end

% Shortcuts to x and v domain ranges (this will go away soon)
Lmin = pde.dimensions{1}.domainMin;
Lmax = pde.dimensions{1}.domainMax;
Vmin = pde.dimensions{2}.domainMin;
Vmax = pde.dimensions{2}.domainMax;

% Level information.
LevX = pde.dimensions{1}.lev;
LevV = pde.dimensions{2}.lev;

% Things to be removed
Deg = pde.dimensions{1}.deg;
Lev = pde.dimensions{1}.lev;
DimX = 1;

% Dimensionality.
Dim = nDims;

% Degree
% pde.params.Deg = Deg;

params = pde.params;

%$
% Set time step.

CFL = 0.1;

dt = Lmax/2^LevX/Vmax/(2*Deg+1)*CFL;

if ~quiet; disp(sprintf('dt = %g', dt )); end

%% Step 1.1. Setup the multi-wavelet transform in 1D (for each dimension).

for d=1:nDims
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde.dimensions{d}.deg,2^pde.dimensions{d}.lev);
end

%% Step 1.2. Apply the mulit-wavelet transform to the initial conditions in each dimension.
% Generate the 1D initial conditions. Input: LevX,LevV,Deg,Lmax,Vmax,pde
% Output: fval (fv and fx)--intial condition f(x,v,t=0)
%         rho--intial condition rho(x,t=0)

if ~quiet; disp('[1.2] Setting up 1D initial conditions'); end
%[fv,fx] = Intial_Con(LevX,LevV,Deg,Lmax,Vmax,pde,FMWT_COMP_x,FMWT_COMP_v);
fx = forwardMWT(pde.dimensions{1}.lev,pde.dimensions{1}.deg,...
    pde.dimensions{1}.domainMin,pde.dimensions{1}.domainMax,...
    pde.dimensions{1}.init_cond_fn,pde.params);
fv = forwardMWT(pde.dimensions{2}.lev,pde.dimensions{2}.deg,...
    pde.dimensions{2}.domainMin,pde.dimensions{2}.domainMax,...
    pde.dimensions{2}.init_cond_fn,pde.params);

%% Step 2. Generate Sparse-Grid (as the Hash + Connectivity tables).

%%% Construct forward and inverse hash tables.
if ~quiet; disp('[2.1] Constructing hash and inverse hash tables'); end
%[HASH,HASHInv,index1D] = HashTable2(Lev,Dim);

[HASH,HASHInv] = HashTable(Lev,Dim,gridType);

nHash = numel(HASHInv);

%%% Construct the connectivity.
if ~quiet; disp('[2.2] Constructing connectivity table'); end
Con2D = Connect2D(Lev,HASH,HASHInv,gridType);

%%% Get the multi-wavelet coefficient representation on the sparse-grid,
%%% i.e., above we transformed each of the initial condition 1D
%%% dependencies into the multi-wavelet basis. Here we combine those N 1D
%%% basis representations into the multi-D (here N=2, so 2D) initial
%%% condition, multi-wavelet representation of the initial condition
%%% specified in PDE.
if ~quiet; disp('[2.3] Calculate 2D initial condition on the sparse-grid'); end
fval = initial_condition_vector(fx,fv,HASHInv,pde);

%% Step 3. Generate time-independent coefficient matrices
% Vlasolv Solver:
%   Operators:
%               vMassV: int_v v*l_i(v)*l_j(v)dv GradV: int_v
%               (l_i(v))'*l_j(v)dv GradX: int_x (m_i(x))'*m_j(x)dx
% Poisson Solver:
%               Operators: DelaX: int_x (m_i(x))''*m_j(x)dx
% Input:
%               LevX, LevV, k, dim, Lmax, Vmax
% Output:
%               2D Matrices--vMassV,GradV,GradX,DeltaX

%% Build the time independent coefficient matricies.
% The original way
if ~quiet; disp('[3.1] Calculate time independent matrix coefficients'); end
[vMassV,GradV,GradX,DeltaX,FluxX,FluxV] = matrix_coeff_TI(LevX,LevV,Deg,Lmin,Lmax,Vmin,Vmax,...
     pde.dimensions{1}.FMWT,pde.dimensions{2}.FMWT);
%%
% The generalized PDE spec way
t = 0;
TD = 0;
pde = getCoeffMats(pde,t,TD);

%%% Generate A_encode / A_data time independent data structures.
if ~quiet; disp('[3.2] Generate A_encode data structure for time independent coefficients'); end
if compression == 3
    % the new matrix construction is as _newCon, only works for 
    % compression= 3
%     A_encode=GlobalMatrixSG(vMassV,GradX,HASHInv,Con2D,Deg);
    A_encode=GlobalMatrixSG_newCon(vMassV,GradX,HASH,Lev,Deg);
else
    % A_data is constructed only once per grid refinement, so can be done
    % on the host side.
    A_data = GlobalMatrixSG_SlowVersion(HASHInv,Con2D,Deg,compression);
end

%% Step 4. Generate time-independent global matrix
% Compute the global matrix for spatial variables "x" by
%
% Poisson Solver: A_Poisson (Hash, Dim_x,k,LevX,DeltaX)
%
% Input: Hash, Dim_x,k,LevX,DeltaX,or nu, eps, CurlCurlX Output: A_Poisson
% Another Idea is to solve Poisson Equation on the finest full grid

if ~quiet; disp('[4] Construct matrix for Poisson solve'); end
if DimX>1
    % Construct DeltaX for DimX
else
    A_Poisson = DeltaX; 
end


%% Step 5. Time Loop
%	Step 5.1 Vlasov Equation
%       Generate time dependent coefficient matrix Generate global matrix
%       A_Vlasov(Hash,coef_mat,Dim) Apply A_Vlasov->f by RK
%	Step 5.2 Poisson Equation
%       Solve Poisson Equation sol_Poisson by A_Poisson(f) Compute
%       E=(sol_Poisson)'


% At time = 0 plotting.
if ~quiet; disp('[5.0] Plotting intial condition'); end


% Construct data for reverse MWT in 2D
[Meval_v,v_node,Meval_x,x_node]=matrix_plot(LevX,LevV,Deg,Lmin,Lmax,Vmin,Vmax,...
    pde.dimensions{1}.FMWT,pde.dimensions{2}.FMWT);
[xx,vv]=meshgrid(x_node,v_node);

% Plot initial condition
if ~quiet
    % Transform from wavelet space to real space
    tmp = Multi_2D(Meval_v,Meval_x,fval,HASHInv,Lev,Deg);
    figure(1000)
    
    f2d0 = reshape(tmp,Deg*2^LevX,Deg*2^LevV)';
    
    ax1 = subplot(1,2,1);
    mesh(xx,vv,f2d0,'FaceColor','interp','EdgeColor','none');
    axis([Lmin Lmax Vmin Vmax])
    %caxis([-range1 +range1]);
    title('df');
    
    ax2 = subplot(1,2,2);
    mesh(xx,vv,f2d0,'FaceColor','interp','EdgeColor','none');
    axis([Lmin Lmax Vmin Vmax])
    %caxis([range2n +range2]);
    title('f');
end

% Write the initial condition to file.
write_fval = 0;
if write_fval; write_fval_to_file(fval,Lev,Deg,0); end

count=1;
plotFreq = 1;
err = 0;

if ~quiet; disp('[7] Advancing time ...'); end
nsteps = max(1,floor( TEND/dt));
for L = 1:nsteps,
    
    tic;
    time(count) = (L-1)*dt;
    timeStr = sprintf('Step %i of %i',L,nsteps);
    
    if ~quiet; disp(timeStr); end
    
    if pde.solvePoisson
        %%% Solve Poisson to get E (from 1-rho=1-int f dv)
        if ~quiet; disp('    [a] Solve poisson to get E'); end
        %[E,u] = PoissonSolve2(LevX,Deg,Lmax,fval,A_Poisson,FMWT_COMP_x,Vmax,index1D);
        [E,u] = PoissonSolve(LevX,Deg,Lmax,fval,A_Poisson,FMWT_COMP_x,Vmax);
    end
    
    if pde.applySpecifiedE
        %%% Apply specified E
        if ~quiet; disp('    [a] Apply specified E'); end
        E = forwardMWT(LevX,Deg,Lmin,Lmax,pde.Ex,pde.params);
        E = E * pde.Et(time(count),params);
    end
    
    Emax = max(abs(Meval_x*E)); % max value on each point for E
    
    %     ax3 = subplot(2,2,3);
    %     plot(x_node,Meval_x*u,'r-o')
    %     title(['time = ',num2str(dt*L)])
    %     ax4 = subplot(2,2,4);
    %     plot(x_node,Meval_x*E,'r-o')
    
    %%% Generate EMassX time dependent coefficient matrix.
    if ~quiet; disp('    [b] Calculate time dependent matrix coeffs'); end
    EMassX = matrix_coeff_TD(LevX,Deg,Lmin,Lmax,E,pde.dimensions{1}.FMWT);
    
    %%
    % Set the dat portion of the EMassX part of E.d_dv term.
    
    pde.terms{2}{1}.dat = E;
    
    %%
    % Now construct the TD coeff_mats.
    
    t = time(count);
    TD = 1;
    pde = getCoeffMats(pde,t,TD);
    
    %% Test new PDE spec based generation of the coeff_matrices
    
    disp( [ 'GradX error : '  num2str(norm(pde.terms{1}{1}.coeff_mat - GradX)/norm(GradX)) ]);
    disp( [ 'vMassV error : ' num2str(norm(pde.terms{1}{2}.coeff_mat - vMassV)/norm(vMassV)) ]);
    disp( [ 'EMassX error : ' num2str(norm(pde.terms{2}{1}.coeff_mat - EMassX)/norm(EMassX)) ]);
    disp( [ 'GradV error : '  num2str(norm(pde.terms{2}{2}.coeff_mat - GradV)/norm(GradV)) ]);
    
    %%% Update A_encode for time-dependent coefficient matricies.
    if ~quiet; disp('    [c] Generate A_encode for time-dependent coeffs'); end
    if compression == 3
    % the new matrix construction is as _newCon, only works for 
    % compression= 3
%         B_encode = GlobalMatrixSG(GradV,EMassX,HASHInv,Con2D,Deg);
        B_encode=GlobalMatrixSG_newCon(GradV,EMassX,HASH,Lev,Deg);
        C_encode=[A_encode B_encode];
    else
        
    end
    
    %%% Advance Vlasov in time with RK3 time stepping method.
    if ~quiet; disp('    [d] RK3 time step'); end
    if compression == 3
        fval = TimeAdvance(C_encode,fval,time(count),dt,compression,Deg,pde,HASHInv);
    else
        
        A_data.GradX     = pde.terms{1}{1}.coeff_mat;
        A_data.vMassV    = pde.terms{1}{2}.coeff_mat;
        A_data.EMassX    = pde.terms{2}{1}.coeff_mat;
        A_data.GradV     = pde.terms{2}{2}.coeff_mat;
        
        %         A_data.vMassV    = vMassV;
        %         A_data.GradX     = GradX;
        %         A_data.GradV     = GradV;
        %         A_data.EMassX    = EMassX;
        
        A_data.FluxX = FluxX;
        A_data.FluxV = FluxV;
        
        % Write the A_data structure components for use in HPC version.
        write_A_data = 0;
        if write_A_data && L==1; write_A_data_to_file(A_data,Lev,Deg); end
        
        fval = TimeAdvance(A_data,fval,time(count),dt,compression,Deg,pde,HASHInv,Vmax,Emax);
        
    end
    
    %%% Write the present fval to file.
    if write_fval; write_fval_to_file(fval,Lev,Deg,L); end
    
    %%% Write data for FK6D test
    
    %     fname = ['tests/vlasov4_time_5_3/fval_',num2str(L,'%3.3i'),'.dat'];
    %     fd = fopen(fname,'w'); % where file.dat is the name you want to save to
    %     fwrite(fd,full(fval),'double'); % where U is the vector/matrix you want to store, double is the typename
    %     fclose(fd);
    
    %%% Plot results
    if mod(L,plotFreq)==0 && ~quiet
        
        figure(1000)
        
        tmp=Multi_2D(Meval_v,Meval_x,fval,HASHInv,Lev,Deg);
        
        f2d = reshape(tmp,Deg*2^LevX,Deg*2^LevV)';
        
        %         ax1 = subplot(1,2,1);
        ax1 = subplot(2,2,1);
        mesh(xx,vv,f2d-f2d0,'FaceColor','interp','EdgeColor','none');
        axis([Lmin Lmax Vmin Vmax])
        %caxis([-range1 +range1]);
        title('df');
        
        %         ax2 = subplot(1,2,2);
        ax2 = subplot(2,2,2);
        mesh(xx,vv,f2d,'FaceColor','interp','EdgeColor','none');
        axis([Lmin Lmax Vmin Vmax])
        %caxis([range2n +range2]);
        title('f');
        
        title(['f @ ', timeStr])
        pause (0.01)
    end
    
    %%% Get the real space solution
    fval_realspace = Multi_2D(Meval_v,Meval_x,fval,HASHInv,Lev,Deg);
    
    %%% Check against known solution
    if pde.checkAnalytic
        
        % Check the wavelet space solution with the analytic solution
        fval_analytic = exact_solution_vector(HASHInv,pde,L*dt);
        err_wavelet = sqrt(mean((fval(:) - fval_analytic(:)).^2));
        disp(['    wavelet space absolute err : ', num2str(err_wavelet)]);
        disp(['    wavelet space relative err : ', num2str(err_wavelet/max(abs(fval_analytic(:)))*100), ' %']);
        
        % Check the real space solution with the analytic solution
        f2d = reshape(fval_realspace,Deg*2^LevX,Deg*2^LevV)';
        f2d_analytic = pde.analytic_solution(xx,vv,L*dt);
        err_real = sqrt(mean((f2d(:) - f2d_analytic(:)).^2));
        disp(['    real space absolute err : ', num2str(err_real)]);
        disp(['    real space relative err : ', num2str(err_real/max(abs(f2d_analytic(:)))*100), ' %']);
        
        err = err_wavelet;
    end
    
    count=count+1;
    t1 = toc;
    disp(['Took ' num2str(t1) ' [s]']);
    
    %%% Save output
    saveOutput = 0;
    if saveOutput
        stat = mkdir('output');
        fName = ['output/f2d-' sprintf('%04.4d',L) '.mat'];
        f2d = reshape(fval_realspace,Deg*2^LevX,Deg*2^LevV)';
        save(fName,'f2d','fval');
    end
    
end

end

