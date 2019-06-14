input_parser = inputParser;

input_parser.KeepUnmatched = true;

default_lev = 3;
default_deg = 2;
default_num_steps = 5;
default_quiet = false;
default_implicit = false;
default_grid_type = 'SG';
valid_grid_types = {'SG','FG'};
check_grid_type = @(x) any(validatestring(x,valid_grid_types));
default_CFL = 0.01;
default_adapt = false;
default_use_oldhash = false;

addRequired(input_parser, 'pde', @isstruct);
addParameter(input_parser,'lev',default_lev, @isnumeric);
addParameter(input_parser,'deg',default_deg, @isnumeric);
addOptional(input_parser,'num_steps',default_num_steps, @isnumeric);
addOptional(input_parser,'quiet',default_quiet,@islogical);
addOptional(input_parser,'implicit',default_implicit, @islogical);
addOptional(input_parser,'grid_type',default_grid_type, check_grid_type);
addOptional(input_parser,'CFL',default_CFL, @isnumeric);
addOptional(input_parser,'adapt',default_adapt, @islogical);
addOptional(input_parser,'use_oldhash',default_use_oldhash, @islogical);

if numel(varargin) == 0 && ~exist('pde','var')
    
    num_parameters = numel(input_parser.Parameters);
    
    disp('ASGarD - Adaptive Sparse Grid Discrization');
    disp('');
    disp('Run with ...');
    disp('');
    disp("asgard(pde_name,'opt_name',opt_val)");
    disp('');
    disp('e.g.,');
    disp('');
    disp("asgard(continuity1,'lev',4,'deg',3,'implicit',true,'CFL',0.1,'adapt',true)");
    disp('');
    disp('Available options are ...');
    
    for p=1:num_parameters
        disp(input_parser.Parameters(p));
    end
    
    error('Run with no arguments, exiting.');
    
end

parse(input_parser,pde,varargin{:})

num_terms       = numel(pde.terms);
num_dimensions  = numel(pde.dimensions);
num_steps       = input_parser.Results.num_steps;

if numel(input_parser.Results.lev) == num_dimensions
    %%
    % Specify lev_vec at runtime to have dimension dependent lev
    for d=1:num_dimensions
        pde.dimensions{d}.lev = input_parser.Results.lev(d);
    end
else
    %%
    % Specify a single lev which all dimensions get
    assert(numel(input_parser.Results.lev)==1);
    for d=1:num_dimensions
        pde.dimensions{d}.lev = input_parser.Results.lev;
    end
end

pde.deg = input_parser.Results.deg;

if isfield(pde,'CFL')
    opts.CFL = pde.CFL;
else
    opts.CFL = input_parser.Results.CFL;
end

opts.quiet = input_parser.Results.quiet;
opts.grid_type = input_parser.Results.grid_type;
opts.implicit = input_parser.Results.implicit;
opts.adapt = input_parser.Results.adapt;
opts.use_oldhash = input_parser.Results.use_oldhash;

if opts.adapt
    opts.use_oldhash = false;
end

if ~isempty(fieldnames(input_parser.Unmatched))
   disp('Extra inputs:')
   disp(input_parser.Unmatched)
   error('Unrecognised input.')
end





