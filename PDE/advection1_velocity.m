function pde = advection1_velocity
%This is a test advection equation to help with the mirror velocity
%equation, which has the general form
%
% df/dt = 1/v^2 d/dv (v^3 f)
%
% Run with
%
% asgard(advection1_velocity,'timestep_method','BE')

pde.CFL = 0.01;

soln = @(v,t) v+t;

dim_v.domainMin = 0.01;
dim_v.domainMax = 10;
dim_v.init_cond_fn = @(v,p,t) soln(v,t);

BCFunc = @(v,t) soln(v,t);

% Domain is (a,b)

% The function is defined for the plane
% x = a and x = b
BCL_fList = { ...
    @(v,p,t) BCFunc(v,t), ... 
    @(t,p) t.*0 + 1
    };

BCR_fList = { ...
    @(v,p,t) BCFunc(v,t), ... % 
    @(t,p) t.*0 + 1
    };

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_v};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE

%% 

% term V1 == 1/v^2 d/dv(v^3  f))
% term V1 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)f(v))   [grad, g2(v) = v^3, BCL= N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^3;

pterm1 = MASS(g1);
pterm2  = GRAD(num_dims,g2,-1,'D','D', BCL_fList, BCR_fList);
termV_s = TERM_1D({pterm1,pterm2});
termV1   = TERM_ND(num_dims,{termV_s});

%%
% Add terms to the pde object

pde.terms = {termV1};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1
source = { ...
    @(v,p,t) 1 - 4*v - 3*t,   ...   % s1v
    @(t,p) t.*0 + 1 ...   % s1t
    };

pde.sources = {source};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...    
    @(v,p,t) soln(v,t), ...
    @(t,p) t.*0 + 1;
    };

%%
% Function to set time step
    function dt=set_dt(pde,CFL)
        
        Lmax = pde.dimensions{1}.domainMax;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end

pde.set_dt = @set_dt;

end