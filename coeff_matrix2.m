%% Construct 1D coefficient matrix
% This routine returns a 2D array representing an operator coefficient
% matrix for a single dimension (1D). Each term in a PDE requires D many coefficient
% matricies. These operators can only use the supported types below.

% function [mat,bcL,bcR] = coeff_matrix2(t,dimension,term_1D,myDim)
function [mat,matD] = coeff_matrix2(pde,t,dim,term_1D)

% Grad
%   \int_T u'v dT = \hat{u}v|_{\partial T} - \int_T uv' dT
% Here we shall choose \hat{u} = AVG(u)+JUMP(u)/2*cval (cval is given)
% need some test

%% Inputs
% t : time
% dimension : an entry in the pde.dimensions array.
% term_1D : one of the dimension entries in an entry in the pde.terms
% array.

%% Supported boundary condition types (BCL and BCR)
% These are chosen in the pde.dimensions options.
% 0 == periodic
% 1 == dirichlet (set value of solution)
% 2 == neumann   (set first derivative of solution)

%% Available coefficient types (coeff_type)
% 1 == G (v .u')   Grad
% 2 == M (v .u )   Mass
% 3 == S (v'.u')   Stiffness

%% Note on global vs local Lax-Friedrichs (LF) flux
% We do not (cannot) use local upwinding or LF because selecting
% either the sign of the flow field or the value of the coefficient C could
% be multivalued within the multi-D solution for a single-D coeff_matrix.
% Examples of value LF:: \hat{f} = \hat{A*u} = {{A*u}}+|A|*(1-C)/2*[[u]]
% Denote LF = (1-C)
%   LF = 0 --> Central Flux
%   LF = 1 --> Upwind Flux
%   LF = Max(df/du) -->Lax-Friedrich Flux

%% TODO ...
% * Choice of flux (may require input C)
% * Other BCs are done, but the RHS (with source) needs more work
% * Picking which term type

%%
% Shortcuts

params = pde.params;

lev = dim.lev;
deg = dim.deg;
xMin = dim.domainMin;
xMax = dim.domainMax;
FMWT = dim.FMWT;
BCL = dim.BCL;
BCR = dim.BCR;
BCL_fList = dim.BCL_fList;
BCR_fList = dim.BCR_fList;

%%
% Shortcuts to term_1d quantities
dat_W = term_1D.dat;
FluxVal = term_1D.LF;
G = term_1D.G;
type = term_1D.type;

%%
% Setup jacobi of variable x and define coeff_mat
N = 2^(lev);
h = (xMax-xMin) / N;
dof_1D = deg * N;

%%
% Set number of quatrature points (should this be order dependant?)
quad_num = 10;

%%
%  Get quadrature points and weights.
%  quad_x(:) is quad_num by 1
%  quad_w(:) is quad_num by 1
[quad_x,quad_w] = lgwt(quad_num,-1,1);

%%
%  Compute the trace values (values at the left and right of each element for all k)
%  p_L(:) is 1 by deg
%  p_R(:) is 1 by deg
p_L = legendre(-1,deg) * 1/sqrt(h);
p_R = legendre(+1,deg) * 1/sqrt(h);

%%
%  Get the basis functions and derivatives for all k
%  p_val(:,:) is quad_num by deg
%  Dp_val(:,:) is quad_num by deg
p_val  = legendre(quad_x,deg)  * 1/sqrt(h);
Dp_val = dlegendre(quad_x,deg) * 1/sqrt(h) * 2/h;


Mass = sparse(dof_1D,dof_1D);
Grad = sparse(dof_1D,dof_1D);
Stif = sparse(dof_1D,dof_1D);
Flux = sparse(dof_1D,dof_1D);
matD = sparse(dof_1D,dof_1D);

%%
% Convert input dat from wavelet (_W) space to realspace (_R)

if isempty(dat_W)
    dat_W = ones(dof_1D,1);
end
dat_R = FMWT' * dat_W;

%% Loop over all elements in this D
%  Here we construct the 1D coeff_mat in realspace, then transform to
%  wavelet space afterwards.
for i=0:N-1
    
    xL = xMin + i*h;
    xR = xL + h;
    
    %%
    % Get index ranges for ...
    
    %%
    %  Current element
    c1 = deg*i+1;
    c2 = deg*(i+1);
    c = c1:c2;
    
    %%
    % Previous element
    p1 = deg*(i-1)+1;
    p2 = deg*i;
    p = p1:p2;
    
    %%
    % Later element
    l1 = deg*(i+1)+1;
    l2 = deg*(i+2);
    l = l1:l2;
    
    %%
    % Map quadrature points from [-1,1] to physical domain of this i element
    
    quad_xi = (((quad_x+1)/2+i)*h+xMin);
    
    %%
    % Build Average (AVG) and Jump (JMP) operators
    
    val_AVG = (1/h) * [-p_L'*p_R/2  -p_L'*p_L/2, ...   % for x1 (left side)
        p_R'*p_R/2   p_R'*p_L/2];      % for x2 (right side)
    
    val_JMP = (1/h) * [ p_L'*p_R    -p_L'*p_L, ...     % for x1 (left side)
        -p_R'*p_R     p_R'*p_L  ];    % for x2 (right side)
    
    %%
    % Combine AVG and JMP to give choice of flux for this operator type
    
    val_FLUX = val_AVG + val_JMP / 2 * FluxVal;
    
    %%
    % Perform volume integral to give deg x deg matrix block
    
    %%
    % Get dat_R at the quadrature points
    
    dat_R_quad = p_val * dat_R(c1:c2);
    
    %%
    % M // mass matrix u . v
    G1 = G(quad_xi,params,t,dat_R_quad);
    val_mass = p_val' * (G1 .* p_val .* quad_w) * h/2;
    
    %%
    % G // grad matrix u . v'
    G1 = G(quad_xi,params,t,dat_R_quad);
    val_grad  = -Dp_val'* (G1 .* p_val .* quad_w) * h/2;
    
    %%
    % S // stiffness matrix u' . v'
    G1 = G(quad_xi,params,t,dat_R_quad);
    val_stif  = Dp_val'* (G1 .* Dp_val .* quad_w) * h/2;
    
    Iu = meshgrid( deg*i+1 : deg*(i+1) );
    
    Mass = Mass + sparse(Iu',Iu,val_mass,dof_1D,dof_1D);
    Grad = Grad + sparse(Iu',Iu,val_grad,dof_1D,dof_1D);
    Stif = Stif + sparse(Iu',Iu,val_stif,dof_1D,dof_1D);
    
    %%
    % Setup numerical flux choice (interior elements only)
    
    if i<N-1 && i>0
        Iu = [ meshgrid(p), meshgrid(c), meshgrid(c), meshgrid(l) ];
        Iv = [ meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
    end
    
    %%
    % Setup boundary conditions
    
    %%
    % If periodic
    
    if BCL == 0 || BCR == 0 %% periodic'
        
        %%
        % Left boundary
        if i==0
            Iu=[meshgrid([deg*(N-1)+1:deg*(N)]),meshgrid(c),meshgrid(c),meshgrid(l)];
            Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
        end
        
        %%
        % Right boundary
        if i==N-1
            Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid([1:deg])];
            Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
        end
        
    end
    
    %%
    % -<funcCoef*{q},p>
    %----------------------------------------------
    % Numerical Flux is defined as
    % Flux = {{f}} + C/2*[[u]]
    %      = ( f_L + f_R )/2 + FunCoef*( u_R - u_L )/2
    % [[v]] = v_R - v_L
    
    FCL = G(xL,params,t,dat_R_quad);
    FCR = G(xR,params,t,dat_R_quad);
    TraVal = [...
        (-p_L)' * FCL * p_R/2 + FluxVal * abs(FCL)/2 * (-p_L)' *   p_R,...
        (-p_L)' * FCL * p_L/2 + FluxVal * abs(FCL)/2 * (-p_L)' * (-p_L),...% xL
        ( p_R)' * FCR * p_R/2 + FluxVal * abs(FCR)/2 * ( p_R)' *   p_R,...
        ( p_R)' * FCR * p_L/2 + FluxVal * abs(FCR)/2 * ( p_R)' * (-p_L),...% xR
        ];
    
    %%
    % If dirichelt
    % u^-_LEFT = g(LEFT)
    % u^+_RIGHT = g(RIGHT)
    
    if BCL == 1 %% left dirichlet
        
        if i==0
            Iu=[meshgrid(c) , meshgrid(c) ,meshgrid(l) ];
            Iv=[meshgrid(c)', meshgrid(c)',meshgrid(c)'];
            
            val_AVG = (1/h) * [...
                -p_L'*(p_L-p_L), ...            % for x1 (left side)
                p_R'*p_R/2   p_R'*p_L/2];      % for x2 (right side)
            
            val_JMP = (1/h) * [...
                -p_L'*(p_L-p_L), ...            % for x1 (left side)
                -p_R'*p_R     p_R'*p_L  ];    % for x2 (right side)
            
            %%
            % Combine AVG and JMP to give choice of flux for this operator type
            
            val_FLUX = val_AVG + val_JMP / 2 * FluxVal;
        end
        
    end
    
    if BCR == 1 %% right dirichlet
        
        if i==N-1
            Iu=[meshgrid(p),meshgrid(c),meshgrid(c)];
            Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)'];
            
            val_AVG = (1/h) * [...
                -p_L'*p_R/2  -p_L'*p_L/2, ...   % for x1 (left side)
                p_R'*(p_R-p_R)];               % for x2 (right side)
            
            val_JMP = (1/h) * [...
                p_L'*p_R    -p_L'*p_L, ...     % for x1 (left side)
                -p_R'*(p_R-p_R)];             % for x2 (right side)
            
            %%
            % Combine AVG and JMP to give choice of flux for this operator type
            
            val_FLUX = val_AVG + val_JMP / 2 * FluxVal;
        end
    end
    
    %%
    % If neumann
    % (gradient u)*n = g
    % by splitting grad u = q by LDG methods, the B.C is changed to
    % q*n = g (=> q = g for 1D variable)
    % only work for derivatives greater than 1
    
    if BCL == 2 %% left neumann
        
        if i==0
            %             Iu=[meshgrid(c) , meshgrid(c) ,meshgrid(l) ];
            %             Iv=[meshgrid(c)', meshgrid(c)',meshgrid(c)'];
            %
            %             val_AVG = (1/h) * [-p_L'*p_L, ...   % for x1 (left side)
            %                 p_R'*p_R/2   p_R'*p_L/2];      % for x2 (right side)
            %
            %             val_JMP = (1/h) * [-p_L'*(p_L-p_L), ...     % for x1 (left side)
            %                 -p_R'*p_R     p_R'*p_L  ];    % for x2 (right side)
            %
            %             %%
            %             % Combine AVG and JMP to give choice of flux for this operator type
            %
            %             val_FLUX = val_AVG + val_JMP / 2 * FluxVal;
                    
            TraVal = [...
                (-p_L)' * (p_L-p_L),...,...
                (-p_L)' * FCL * p_L,...% xL
                ( p_R)' * FCR * p_R/2 + FluxVal * abs(FCR)/2 * ( p_R)' *   p_R,...
                ( p_R)' * FCR * p_L/2 + FluxVal * abs(FCR)/2 * ( p_R)' * (-p_L),...% xR
                ];
            
        end
    end
    
    if BCR == 2 %% right neumann
        
        if i==N-1
            %             Iu=[meshgrid(p),meshgrid(c),meshgrid(c)];
            %             Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)'];
            %
            %             val_AVG = (1/h) * [-p_L'*p_R/2  -p_L'*p_L/2, ...   % for x1 (left side)
            %                 p_R'*p_R];      % for x2 (right side)
            %
            %             val_JMP = (1/h) * [ p_L'*p_R    -p_L'*p_L, ...     % for x1 (left side)
            %                 -p_R'*(p_R-p_R)];    % for x2 (right side)
            %
            %             %%
            %             % Combine AVG and JMP to give choice of flux for this operator type
            %
            %             val_FLUX = val_AVG + val_JMP / 2 * FluxVal;
            %
            
            TraVal = [...
                (-p_L)' * FCL * p_R/2 + FluxVal * abs(FCL)/2 * (-p_L)' *   p_R,...
                (-p_L)' * FCL * p_L/2 + FluxVal * abs(FCL)/2 * (-p_L)' * (-p_L),...% xL
                ( p_R)' * FCR * p_R,...
                ( p_R)' * (p_R-p_R),...% xR
                ];
        end
    end
    
    
    %%
    % Adding trace value to matrix
    
    RowInd = [c'*ones(1,deg) c'*ones(1,deg) c'*ones(1,deg) c'*ones(1,deg)];
    ColInd = [ones(deg,1)*(c-deg),ones(deg,1)*c,ones(deg,1)*c,ones(deg,1)*(c+deg)];
    
    if i == 0
        Iu = RowInd(:,deg+1:end);
        Iv = ColInd(:,deg+1:end);
        Val = TraVal(:,deg+1:end);
    elseif i == N - 1
        Iu = RowInd(:,1:3*deg);
        Iv = ColInd(:,1:3*deg);
        Val = TraVal(:,1:3*deg);
    else
        Iu = RowInd;
        Iv = ColInd;
        Val = TraVal;
    end
    
    %     Grad = Grad - sparse(Iv,Iu,val_FLUX,dof_1D,dof_1D);
    Grad = Grad + sparse(Iu,Iv,Val,dof_1D,dof_1D);
    
end


%% Transform coeff_mat to wavelet space
Mass = FMWT * Mass * FMWT';
Grad = FMWT * Grad * FMWT';

if type == 3
    
    % Use LDG method, i.e., split into two first order equations, then
    % recombine
    
    %%
    % Get a grad operator with downwind flux and Neumann BCs
    
    termA = term_1D;
    dimA = dim;
    
    termA.type = 1;      % Grad Operator
    termA.LF = -1;       % Downwind Flux
    termA.G = @(x,p,t,dat) x*0+1;
    dimA.BCL = 2; % Neumann
    dimA.BCR = 2; % Neumann
    matD = coeff_matrix2(pde,t,dimA,termA);
    
    %%
    % Get a grad operator with upwind flux and Dirichlet BCs
    
    termB = term_1D;
    dimB = dim;
    
    termB.type = 1;      % Grad Operator
    termB.LF = 1;       % Upwind Flux
    %     termB.G = @(x,p,t,dat) x*0+1;
    dimB.BCL = 1; % Dirichlet
    dimB.BCR = 1; % Dirichlet
    matU = coeff_matrix2(pde,t,dimB,termB);
    
    %%
    % Combine back into second order operator
    Stif = matD*matU;
    
end

if type == 1
    mat = Grad;
end
if type == 2
    mat = Mass;
end
if type == 3
    mat = Stif;
end

end