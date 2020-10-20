function pde = seperable_adapt_test_2d(opts)
% 2D test case using continuity equation, i.e.,
%
% df/dt == -df/dx
%
% Run with
%
% asgard(@seperable_adapt_test_2d,'timestep_method','BE','adapt',true,'many_solution_capable',true,'adapt_initial_condition',true,'adapt_threshold',1e-1)

w = 1.0;
C = 1.0;

    function ret = offset(t)
        ret = sin(t);
    end
    function ret = expf(x,t)
        ret = exp( -( x - offset(t) ).^2 ./ w );
    end


%% Define the analytic solution (optional).
% This requires nDims+time function handles.

solution1 = { ...
    @(x,p,t)  x .* expf(x,t), ...
    @(y,p,t)  y.*0 + 1, ...
    @(t)      1 ...
    };

solution2 = { ...
    @(x,p,t)  expf(x,t), ...
    @(y,p,t)  -C.*sin(2*pi*y), ...
    @(t)      1 ...
    };

solutions = {solution1,solution2};

%%
% Add initial conditions

initial_conditions = {solution1,solution2};

%% Define the dimensions
%
% Here we setup a 2D problem (x,y)

dim_x = DIMENSION(-10,+10);
dim_y = DIMENSION(-1,+1);

dimensions = {dim_x,dim_y};
num_dims = numel(dimensions);

%% Define the terms of the PDE
%
% Here we have 2 terms, having only nDims=2 (x,y) operators.

%%
% -df/dx which is 
%
% d/dx g1(x) f(x,y)          [grad,g1(x)=-1, BCL=P, BCR=P]

g1 = @(x,p,t,dat) x*0-1;
pterm1  = GRAD(num_dims,g1,0,'D','D');
term1_x = TERM_1D({pterm1});
term1   = TERM_ND(num_dims,{term1_x,[]});

terms = {term1};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

source1_t1p1 = { ...
    @(x,p,t) 2./w.*expf(x,t).*x.^2*cos(t), ... 
    @(y,p,t) y.*0+1, ...
    @(t)     1 ...
    };

source1_t1p2 = { ...
    @(x,p,t)  -2/w.*expf(x,t).*x*cos(t)*sin(t), ...
    @(y,p,t)  y.*0+1, ...
    @(t)      1 ...
    };

source1_t2p1 = { ...
    @(x,p,t)  expf(x,t), ...
    @(y,p,t)  y.*0+1, ...
    @(t)      1 ...
    };

source1_t2p2 = { ...
    @(x,p,t)  -2/w.*x.^2.*expf(x,t), ...
    @(y,p,t)  y.*0+1, ...
    @(t)      1 ...
    };

source1_t2p3 = { ...
    @(x,p,t)  2/w.*x.*expf(x,t).*sin(t), ...
    @(y,p,t)  y.*0+1, ...
    @(t)      1 ...
    };

source2_t1p1 = { ...
    @(x,p,t)  -2*C/w.*expf(x,t).*x.*cos(t), ...
    @(y,p,t)  sin(2*pi*y), ...
    @(t)      1 ...
    };

source2_t1p2 = { ...
    @(x,p,t)  2*C/w.*expf(x,t).*cos(t).*sin(t), ...
    @(y,p,t)  sin(2*pi*y), ...
    @(t)      1 ...
    };

source2_t2p1 = { ...
    @(x,p,t)  2*C/w.*expf(x,t).*x, ...
    @(y,p,t)  sin(2*pi*y), ...
    @(t)      1 ...
    };

source2_t2p2 = { ...
    @(x,p,t)  -2*C/w.*expf(x,t).*sin(t), ...
    @(y,p,t)  sin(2*pi*y), ...
    @(t)      1 ...
    };

sources = {source1_t1p1,source1_t1p2,...
    source1_t2p1,source1_t2p2,source1_t2p3...
    source2_t1p1,source2_t1p2,...
    source2_t2p1,source2_t2p2};

%% Define function to set dt

    function dt=set_dt(pde,CFL)       
        Lmax = pde.dimensions{1}.max;
        Lmin = pde.dimensions{1}.min;
        LevX = pde.dimensions{1}.lev;
        dt = (Lmax-Lmin)/2^LevX*CFL;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end


