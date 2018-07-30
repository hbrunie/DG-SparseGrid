function [GradMat,GradGradMat] = Matrix_TI(Lev,Deg,Lmax,pde)
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
dofs = 3*dof_1D_x^3;
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

% Matrix for (u,v')
val = Dp_val'*(quad_w*ones(1,2).*p_val);
Ac = repmat({val},nx,1);
G = blkdiag(Ac{:});

% Matrix for (u',v)
val = p_val'*(quad_w*ones(1,2).*Dp_val);
Ac = repmat({val},nx,1);
G2 = blkdiag(Ac{:});


% Matrix for (u,v) = eye


%****************************************
% Term for <{u_h},[v]>ds
%****************************************
Amd  = -p_1'*p_1/2+p_2'*p_2/2;
Asub = -p_1'*p_2/2;
Asup =  p_2'*p_1/2;
K = 1/hx*blktridiag([Amd],[Asub],[Asup],nx);

%****************************************
% Term for <[u_h],{v}>ds
%****************************************
Amd  = p_1'*(-p_1)/2+p_2'*p_2/2;
Asub =  p_1'*p_2/2;
Asup =  p_2'*(-p_1)/2;
H = 1/hx*blktridiag([Amd],[Asub],[Asup],nx);

%****************************************
% Term for <{u_h'},[v]>ds
%****************************************
Amd = -p_1'*Dp_1/2+p_2'*Dp_2/2;
Asub =-p_1'*Dp_2/2;
Asup = p_2'*Dp_1/2;
L = blktridiag([Amd]/hx,[Asub]/hx,[Asup]/hx,nx);

%****************************************
% Term for <[u_h],{v'}>ds
%****************************************
Amd = Dp_1'*(-p_1)/2+Dp_2'*p_2/2;
Asub = Dp_1'*p_2/2;
Asup = Dp_2'*(-p_1)/2;
J = blktridiag([Amd]/hx,[Asub]/hx,[Asup]/hx,nx);

%****************************************
% Term for <[u_h],[v]>ds
%****************************************
Amd  =  p_1'*(p_1)+p_2'*p_2;
Asub = -p_1'*p_2/2;
Asup =  p_2'*(-p_1)/2;
Q = blktridiag([Amd]/hx,[Asub]/hx,[Asup]/hx,nx);

II = speye(dof_1D_x);
T1=[...
    kron(kron(II,GG),II)+kron(kron(II,II),GG),-kron(kron(G,G'),II),-kron(kron(G,II),G');...
    -kron(kron(G',G),II),kron(kron(GG,II),II)+kron(kron(II,II),GG),-kron(kron(II,G),G');...
    -kron(kron(G',II),G),-kron(kron(II,G'),G),kron(kron(GG,II),II)+kron(kron(II,GG),II);...
];
T2=[...
    kron(kron(II,J),II)+kron(kron(II,II),J),-kron(kron(H,G'),II),-kron(kron(H,II),G');...
    -kron(kron(G',H),II),kron(kron(J,II),II)+kron(kron(II,II),J),-kron(kron(II,H),G');...
    -kron(kron(G',II),H),-kron(kron(II,G'),H),kron(kron(J,II),II)+kron(kron(II,J),II);...
];
T3=[...
    kron(kron(II,L),II)+kron(kron(II,II),L),-kron(kron(G,K),II),-kron(kron(G,II),K);...
    -kron(kron(K,G),II),kron(kron(L,II),II)+kron(kron(II,II),L),-kron(kron(II,G),K);...
    -kron(kron(K,II),G),-kron(kron(II,K),G),kron(kron(L,II),II)+kron(kron(II,L),II);...
];
T4=[kron(kron(II,Q),II)+kron(kron(II,II),Q),zeros(dof_1D_x^3),zeros(dof_1D_x^3);...
    zeros(dof_1D_x^3),kron(kron(Q,II),II)+kron(kron(II,II),Q),zeros(dof_1D_x^3);...
    zeros(dof_1D_x^3),zeros(dof_1D_x^3),kron(kron(Q,II),II)+kron(kron(II,Q),II)];

% generate A_encode
[IndexI,IndexJ]=meshgrid([1:dof_1D_x^3]);
% A11
count = 1;
A_encode{count}.IndexI = IndexI';
A_encode{count}.IndexJ = IndexJ';

A_encode{count}.A1=II;
A_encode{count}.A2=G-J-L+alpha*Q;
A_encode{count}.A3=II;

count = count+1;
A_encode{count}.IndexI = IndexI';
A_encode{count}.IndexJ = IndexJ';

A_encode{count}.A1=II;
A_encode{count}.A2=II;
A_encode{count}.A3=G-J-L+alpha*Q;

% A12
count = count+1;
A_encode{count}.IndexI = IndexI';
A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';

A_encode{count}.A1=-G;
A_encode{count}.A2=G';
A_encode{count}.A3=II;

count = count+1;
A_encode{count}.IndexI = IndexI';
A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';

A_encode{count}.A1=H;
A_encode{count}.A2=G';
A_encode{count}.A3=II;

count = count+1;
A_encode{count}.IndexI = IndexI';
A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';

A_encode{count}.A1=G;
A_encode{count}.A2=H';
A_encode{count}.A3=II;

% A13
count = count+1;
A_encode{count}.IndexI = IndexI';
A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';

A_encode{count}.A1=-G;
A_encode{count}.A2=II;
A_encode{count}.A3=G';

count = count+1;
A_encode{count}.IndexI = IndexI';
A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';

A_encode{count}.A1=H;
A_encode{count}.A2=II;
A_encode{count}.A3=G';

count = count+1;
A_encode{count}.IndexI = IndexI';
A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';

A_encode{count}.A1=G;
A_encode{count}.A2=II;
A_encode{count}.A3=H';

% A21
count = count+1;
A_encode{count}.IndexI = dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = IndexJ';

A_encode{count}.A1=-G';
A_encode{count}.A2=G;
A_encode{count}.A3=II;

count = count+1;
A_encode{count}.IndexI = dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = IndexJ';

A_encode{count}.A1=G';
A_encode{count}.A2=H;
A_encode{count}.A3=II;

count = count+1;
A_encode{count}.IndexI = dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = IndexJ';

A_encode{count}.A1=H';
A_encode{count}.A2=G;
A_encode{count}.A3=II;

% A22
count = count+1;
A_encode{count}.IndexI = dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';

A_encode{count}.A1=G-J-L+alpha*Q;
A_encode{count}.A2=II;
A_encode{count}.A3=II;

count = count+1;
A_encode{count}.IndexI = dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';

A_encode{count}.A1=II;
A_encode{count}.A2=II;
A_encode{count}.A3=G-J-L+alpha*Q;

% A23
count = count+1;
A_encode{count}.IndexI = dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';

A_encode{count}.A1=-II;
A_encode{count}.A2=G;
A_encode{count}.A3=G';

count = count+1;
A_encode{count}.IndexI = dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';

A_encode{count}.A1=II;
A_encode{count}.A2=H;
A_encode{count}.A3=G';

count = count+1;
A_encode{count}.IndexI = dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';

A_encode{count}.A1=II;
A_encode{count}.A2=G;
A_encode{count}.A3=H';

% A31
count = count+1;
A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = IndexJ';

A_encode{count}.A1=-G';
A_encode{count}.A2=II;
A_encode{count}.A3=G;

count = count+1;
A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = IndexJ';

A_encode{count}.A1=G';
A_encode{count}.A2=II;
A_encode{count}.A3=H;

count = count+1;
A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = IndexJ';

A_encode{count}.A1=H';
A_encode{count}.A2=II;
A_encode{count}.A3=G;

% A32
count = count+1;
A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';

A_encode{count}.A1=-II;
A_encode{count}.A2=G';
A_encode{count}.A3=G;

count = count+1;
A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';

A_encode{count}.A1=II;
A_encode{count}.A2=G';
A_encode{count}.A3=H;

count = count+1;
A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = dof_1D_x^3+IndexJ';

A_encode{count}.A1=II;
A_encode{count}.A2=H';
A_encode{count}.A3=G;

% A33
count = count+1;
A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';

A_encode{count}.A1=G-J-L+alpha*Q;
A_encode{count}.A2=II;
A_encode{count}.A3=II;

count = count+1;
A_encode{count}.IndexI = 2*dof_1D_x^3+IndexI';
A_encode{count}.IndexJ = 2*dof_1D_x^3+IndexJ';

A_encode{count}.A1=II;
A_encode{count}.A2=G-J-L+alpha*Q;
A_encode{count}.A3=II;

% size(T1)
% dofs
Mat = T1-T2-T3+alpha*T4-pde.w2*speye(dofs,dofs);
figure;subplot(2,2,1);spy(T1);subplot(2,2,2);spy(T2);subplot(2,2,3);spy(T3);subplot(2,2,4);spy(T4);
figure;spy(Mat);

% handling boundary condition

% condest(Mat)

Lend = Lmax;Lstart = 0;
fx = zeros(quad_num,3);
fy = zeros(quad_num,3);
fz = zeros(quad_num,3);
for i=0:2^Lev-1
    hx = (Lend-Lstart)/nx;
    xi_x = hx*(quad_x/2+1/2+i)+Lstart;
    fx = pde.rhs(xi_x,1/2+xi_x-xi_x,1/2+xi_x-xi_x);
    fy = pde.rhs(1/2+xi_x-xi_x,xi_x,1/2+xi_x-xi_x);
    fz = pde.rhs(1/2+xi_x-xi_x,1/2+xi_x-xi_x,xi_x);
    
    index = Deg*i+1:Deg*(i+1);
    f1x(index) = p_val'*(quad_w.*fx(:,1))*hx*sqrt(1/hx)/2;
    f1y(index) = p_val'*(quad_w.*fy(:,1))*hx*sqrt(1/hx)/2;
    f1z(index) = p_val'*(quad_w.*fz(:,1))*hx*sqrt(1/hx)/2;
%     f1 = kron(kron(f1x,f1y),f1z);


    f2x(index) = p_val'*(quad_w.*fx(:,2))*hx*sqrt(1/hx)/2;
    f2y(index) = p_val'*(quad_w.*fy(:,2))*hx*sqrt(1/hx)/2;
    f2z(index) = p_val'*(quad_w.*fz(:,2))*hx*sqrt(1/hx)/2;
%     f2 = kron(kron(f2x,f2y),f2z);

    f3x(index) = p_val'*(quad_w.*fx(:,3))*hx*sqrt(1/hx)/2;
    f3y(index) = p_val'*(quad_w.*fy(:,3))*hx*sqrt(1/hx)/2;
    f3z(index) = p_val'*(quad_w.*fz(:,3))*hx*sqrt(1/hx)/2;
%     f3 = kron(kron(f3x,f3y),f3z);

end

ff =[kron(kron(f1x,f1y),f1z),kron(kron(f2x,f2y),f2z),kron(kron(f3x,f3y),f3z)]';
% size(ff)
% size(Mat)
sol = Mat\ff;
% full([sol,ff/(2*pi^2-1)])
figure;plot(sol-ff/(2*pi^2-1))
max(abs(sol-ff/(2*pi^2-1)))/norm(ff/(2*pi^2-1))
% compute the RHS term

return


GradX = FMWT_COMP_x*GradX*FMWT_COMP_x';

end