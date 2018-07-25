function [GradMat,GradGradMat] = Matrix_TI(Lev,Deg,Lmax,FMWT_COMP_x)
%========================================================
% Construct the matrix for curl operator on [0,Lmax]
% Input: 
%   Lev denotes the Level for the mesh
%   Deg denotes the degree for polynomial
%   Lmax denotes the domain for x variable
%   FMWT_COMP_x denotes the converting relationship
% Output:
%   GradMat:  \int u'v dx
%   GradGradMat: \int u'v' dx
%========================================================
%--DG parameters
quad_num=10;
%---------------

alpha = 10;


% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre(1,Deg);

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);
Dp_1 = dlegendre(-1,Deg);
Dp_2 = dlegendre(1,Deg);


%---------------------------
% Define Matrices
%---------------------------
nx = 2^(Lev);
hx = Lmax/nx;
dof_1D_x = Deg*nx;
GradMat = sparse(dof_1D_x,dof_1D_x);
GradGradMat = sparse(dof_1D_x,dof_1D_x);

%******************************************************
% generate 1D matrix for DG (du/dx,v) by weak form
% Here we only consider the central flux in the scheme
% -(u,dv/dx)+<{u},[v]>
%******************************************************
%======================================
% Matrices related to x variable GradX
%======================================
% compute (u',v)+1/2*u^{-}[v^{+}-v^{-}]
% generate 1D matrix for DG

% Matrix for (u',v')
val = Dp_val'*(quad_w*ones(1,2).*Dp_val)*1/hx;
Ac = repmat({val},nx,1);
GG = blkdiag(Ac{:});

% Matrix for (u',v)
val = Dp_val'*(quad_w*ones(1,2).*p_val);
Ac = repmat({val},nx,1);
G = blkdiag(Ac{:});

% Matrix for (u,v) = eye


%****************************************
% Term for <{u_h},[v]>ds
% Numerical Flux is taken as central flux
%****************************************
Amd  = -p_1'*p_1/2+p_2'*p_2/2;
Asub = -p_1'*p_2/2;
Asup =  p_2'*p_1/2;
T1 = blktridiag([Amd]/hx,[Asub]/hx,[Asup]/hx,nx);

%****************************************
% Term for <[u_h],[v]>ds
% Numerical Flux is taken as central flux
%****************************************
Amd =  p_1'*p_1+p_2'*p_2;
Asub = -p_1'*p_2;
Asup = -p_2'*p_1;
T2 = blktridiag([Amd]/hx,[Asub]/hx,[Asup]/hx,nx)*alpha;

%****************************************
% Term for <{u_h'},[v]>ds
% Numerical Flux is taken as central flux
%****************************************
Amd = -p_1'*Dp_1/2+p_2'*Dp_2/2;
Asub =-p_1'*Dp_2/2;
Asup = p_2'*Dp_1/2;
T3 = blktridiag([Amd]/hx,[Asub]/hx,[Asup]/hx,nx);

%****************************************
% Term for <[u_h],{v'}>ds
% Numerical Flux is taken as central flux
%****************************************
Amd = -Dp_1'*p_1/2+Dp_2'*p_2/2;
Asub =  Dp_1'*p_2/2;
Asup = -Dp_2'*p_1/2;
T4 = blktridiag([Amd]/hx,[Asub]/hx,[Asup]/hx,nx);

S = GG - T3 -T4+T2;
H = -G'+T1;


figure;subplot(2,2,1);spy(T1);subplot(2,2,2);spy(T2);subplot(2,2,3);spy(T3);subplot(2,2,4);spy(T4);
figure;subplot(1,2,1);spy(S);subplot(1,2,2);spy(H);



GradX = FMWT_COMP_x*GradX*FMWT_COMP_x';

end