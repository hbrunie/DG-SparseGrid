function pde = advection1(opts)
% 1D test case using continuity equation, i.e., 
%
% df/dt == -2*df/dx - 2*sin(x)
%
% Run with
%
% explicit
% asgard(@advection1)
% asgard(@advection1,'lev',4,'deg',3)
%
% implicit
% asgard(@advection1,'timestep_method','CN')
% asgard(@advection1,'timestep_method','CN','CFL',0.01)


%% Define dimensions

dim_x = DIMENSION(0,pi);
dim_x.init_cond_fn = @(x,p,t) cos(x);

dimensions = {dim_x};

num_dims = numel(dimensions);


%% Define boundary conditions

% Dirichlet on the right, f(0) = 1

BCFunc_Left  = @(x) x.*0 +1;
BCFunc_Right = @(x) x.*0 -1;

BCL_fList = { ...
    @(x,p,t) BCFunc_Left(x), ... % replace x by b
    @(t,p) t.*0 + 1
    };

BCR_fList = { ...
    @(x,p,t) BCFunc_Right(x), ... % replace x by b
    @(t,p) t.*0 + 1
    };


%% Define PDE terms
% Here we have 1 term1, having only nDims=1 (x) operators.
 
% -2*df/dx
g1 = @(x,p,t,dat) x.*0-2;
pterm1 = GRAD(num_dims,g1,-1,'D','D', BCL_fList, BCR_fList);

term1_x = TERM_1D({pterm1});
term1   = TERM_ND(num_dims,{term1_x});

terms = {term1};


%% Define paramaters
%  These might be used within the various functions.

params.parameter1 = 0;
params.parameter2 = 1;


%% Define sources

% Source 1
s1x = @(x,p,t) -2.*sin(x);
s1t = @(t,p) t.*0 + 1;
source1 = {s1x,s1t};

sources = {source1};


%% Define the analytic solution (optional).
% This requires nDims+time function handles.

a_x = @(x,p,t) cos(x);
a_t = @(t,p) t.*0 + 1;

analytic_solutions_1D = {a_x,a_t};

%% Define function to calculate time step

    function dt=set_dt(pde,CFL)
        
        dim = pde.dimensions{1};
        lev = dim.lev;
        xMax = dim.max;
        xMin = dim.min;
        xRange = xMax-xMin;
        dx = xRange/(2^lev);
        dt = CFL*dx;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,analytic_solutions_1D);

end
