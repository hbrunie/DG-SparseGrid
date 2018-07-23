function Matrix=matrix_coeff_TI_v3(Lev,k,Lmax,Vmax)
%=============================================================
% Algorithm 3:: Time-independent Matrices
% Generate time-independent coefficient matrices
% Vlasolv Solver:
%   Operators:  vMassV: int_v v*l_i(v)*l_j(v)dv
%               GradV: int_v (l_i(v))'*l_j(v)dv
%               GradX: int_x (m_i(x))'*m_j(x)dx
%               NGradX: set numerical flux as 
% 
% Input: Lev, k, dim, Lmax, Vmax
% Choose Cval for upwinding or alternating flux
% P.S. This is the full-grid version
%=============================================================
Matrix = struct();

load(['two_scale_rel_',num2str(k),'.mat'])
n=Lev;

H0(find(abs(H0)<1e-5))=0;
G0(find(abs(G0)<1e-5))=0;

H1 = zeros(k);
G1 = zeros(k);

for j_x = 1:k
    for j_y = 1:k
        H1(j_x,j_y) = ((-1)^(j_x+j_y-2)  )*H0(j_x,j_y);
        G1(j_x,j_y) = ((-1)^(k+j_x+j_y-2))*G0(j_x,j_y);
    end
end

FMWT = zeros(k*2^Lev);
iFMWT = zeros(k*2^Lev);

Np=2^Lev;
for j=1:2^Lev/2
    % The reverse order from Lin
    FMWT(k*(j-1)+1:k*j,2*k*(j-1)+1:2*k*j)=[H0 H1];
    FMWT(k*(j+Np/2-1)+1:k*(j+Np/2),2*k*(j-1)+1:2*k*j) = [G0 G1];
end
iFMWT=FMWT';

sp = [];
FMWT_COMP = eye(k*Np);
for j=1:n
    cFMWT = FMWT;
    % Reverse the index in matrix from Lin
    if j>1
        cFMWT = zeros(k*Np);
        cn = 2^(n-j+1)*k;
        cnr=Np*k-cn;
        cFMWT(cn+1:k*Np,cn+1:k*Np)=eye(Np*k-cn);
        cFMWT(1:cn/2,1:cn)=FMWT(1:cn/2,1:cn);
        cFMWT(cn/2+1:cn,1:cn)=FMWT(k*Np/2+1:k*Np/2+cn/2,1:cn);
    end

    FMWT_COMP = cFMWT*FMWT_COMP;
end



%--Quadrature
quad_num=10;
%---------------

% compute the trace values
p_1 = legendre(-1,k);
p_2 = legendre(1,k);

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,k);
Dp_val = dlegendre(quad_x,k);

%---------------------------
% Jacobi of variable x and v
% Define Matrices
%---------------------------
nx=2^(Lev);hx=Lmax/nx;
Jacobi_x=hx;
dof_1D_x=k*nx;
GradX=sparse(dof_1D_x,dof_1D_x);
DeltaX=sparse(2*dof_1D_x,2*dof_1D_x);


nv=2^(Lev);hv=2*Vmax/nv;
Jacobi_v=hv;
dof_1D_v=k*nv;
vMassV=sparse(dof_1D_v,dof_1D_v);
GradV=sparse(dof_1D_v,dof_1D_v);

%***************************
% Term for \int_K(u_h)(v')dK
%***************************
val = 1/hx*[Dp_val'*(quad_w.*p_val)];
Ac = repmat({val},nx,1);
GradX = blkdiag(Ac{:});
full(GradX)

GradX  = FMWT_COMP*GradX*FMWT_COMP';
full(GradX)

H = [H0 H1];
G = [G0 G1];

A{1,1} = [Dp_val'*(quad_w.*p_val)]*2^0;

% A{1,2} = H*blkdiag(A{1,1}*2,A{1,1}*2)*G';
% A{2,1} = G*blkdiag(A{1,1}*2,A{1,1}*2)*H';
% A{2,2} = G*blkdiag(A{1,1}*2,A{1,1}*2)*G';

A{1,2} = H*kron(eye(2),A{1,1})*G'*2;
A{2,1} = G*kron(eye(2),A{1,1})*H'*2;
A{2,2} = G*kron(eye(2),A{1,1})*G'*2;

A{1,3} = H*(kron(A{1,2},eye(2)))*2;
A{2,3} = G*(kron(A{1,2},eye(2)))*2;
A{3,3} = (kron(eye(2),A{2,2}))*2;
A{3,1} = kron(eye(2),A{2,1})*H'*2;
A{3,2} = kron(eye(2),A{2,1})*G'*2;

A{1,4} = H*(kron(A{1,3},eye(2)))*2;
A{2,4} = G*(kron(A{1,3},eye(2)))*2;
A{3,4} = (kron(eye(2),A{2,3}))*2;
A{4,4} =  (kron(eye(2),A{3,3}))*2;
A{4,1} = kron(eye(2),A{3,1})*H'*2;
A{4,2} = kron(eye(2),A{3,1})*G'*2;
A{4,3} = kron(eye(2),A{3,2})*2;

A{1,5} = H*(kron(A{1,4},eye(2)))*2;
A{2,5} = G*(kron(A{1,4},eye(2)))*2;
A{3,5} = (kron(eye(2),A{2,4}))*2;
A{4,5} =  (kron(eye(2),A{3,4}))*2;
A{5,5} =  (kron(eye(2),A{4,4}))*2;
A{5,1} = kron(eye(2),A{4,1})*H'*2;
A{5,2} = kron(eye(2),A{4,1})*G'*2;
A{5,3} = kron(eye(2),A{4,2})*2;
A{5,4} = kron(eye(2),A{4,3})*2;

B = [A{1,1} A{1,2} A{1,3} A{1,4} A{1,5};...
    A{2,1} A{2,2} A{2,3} A{2,4} A{2,5};...
    A{3,1} A{3,2} A{3,3} A{3,4} A{3,5};...
    A{4,1} A{4,2} A{4,3} A{4,4} A{4,5};...
    A{5,1} A{5,2} A{5,3} A{5,4} A{5,5};...
    ]

plot(B-GradX)
return
%****************************************
% Term for \int_{\partial K} {{u_h}}vds
% Numerical Flux is taken as central flux
%****************************************
Amd  = -p_1'*p_1/2/hx+p_2'*p_2/2/hx;
Asub = -p_1'*p_2/2/hx;
Asup = zeros(k,k)+p_2'*p_1/2/hx;
GradXFluxC = -blktridiag([Amd],[Asub],[Asup],nx);
% Adding Periodic Boundary Conditions
IndexEnd = dof_1D_x-k+1:dof_1D_x;
IndexSta = 1:k;
Iu = repmat(IndexEnd,k,1);
Iv = repmat(IndexSta,k,1);
GradXFluxC = GradXFluxC...
    +sparse([Iv',Iu'],[Iu,Iv],-[Asub,Asup],dof_1D_x,dof_1D_x);

%**************************************
% Term for \int_{\partial K} [[u_h]]vds
% Numerical Flux is taken as jump
%**************************************
Amd  = -p_1'*p_1/hx+p_2'*(-p_2)/hx;
Asub = -p_1'*(-p_2)/hx;
Asup = zeros(k,k)+p_2'*p_1/hx;
GradXFluxJ = -blktridiag([Amd],[Asub],[Asup],nx);

% Adding Periodic Boundary Conditions
IndexEnd = dof_1D_x-k+1:dof_1D_x;
IndexSta = 1:k;
Iu = repmat(IndexEnd,k,1);
Iv = repmat(IndexSta,k,1);
GradXFluxJ = GradXFluxJ...
    +sparse([Iv',Iu'],[Iu,Iv],-[Asub,Asup],dof_1D_x,dof_1D_x);


% PGradX denote the flux is taken as \hat{u}=u+
% NGradX denote the flux is taken as \hat{u}=u-
% GradX denote the flux is taken as \hat{u}={u}
PGradX = GradX+GradXFluxC+GradXFluxJ/2;
NGradX = GradX+GradXFluxC-GradXFluxJ/2;
 GradX = GradX+GradXFluxC;

%======================================
% Matrices related to v variable
% vMassV and GradV
%======================================
GradV = (GradX)/2*Lmax/Vmax;

% (vf(v),w(v))_Kv
for Lv=0:nv-1
    xi_v=(( (quad_x+1)/2+Lv)*hv-Vmax); % mapping from [-1,1] to physical domain
    %---------------------------------------------
    % Matrix
    %---------------------------------------------
    % value of local matrix
    val_loc=p_val'*(p_val.*xi_v.*quad_w)*Jacobi_v/2/hv;
    Iu=meshgrid(k*Lv+1:k*(Lv+1));
    vMassV=vMassV+sparse(Iu',Iu,val_loc,dof_1D_v,dof_1D_v);
    
end


%***************************************
% Following is for Multiwavelet DG Basis
% Simple way to convert??
%***************************************
% Transfer to multi-DG bases
vMassV = FMWT_COMP_v*vMassV*FMWT_COMP_v';
GradX  = FMWT_COMP*GradX*FMWT_COMP';
NGradX = FMWT_COMP_x*NGradX*FMWT_COMP_x';
PGradX = FMWT_COMP_x*PGradX*FMWT_COMP_x';

GradV = FMWT_COMP_v*GradV*FMWT_COMP_v';

% % DeltaX = blkdiag(FMWT_COMP_x,FMWT_COMP_x)*...
% %                 DeltaX*...
% %          blkdiag(FMWT_COMP_x',FMWT_COMP_x');

Matrix = struct('GradX',GradX,'NGradX',NGradX,'PGradX',PGradX,...
                'GradV',GradV,'vMassV',vMassV);
