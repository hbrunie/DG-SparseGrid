% main code for LDG 
% of Full Pitch Angle Dynamics ::
% df/dt = - E*d/dx[(1-x^2)f]+C*d/dx[(1-x^2)df/dx]-R*d/dx[x(1-x^2)f]
% with f(t=0)=f0(x), f(x=+,-1)=0.

clear all
close all
% clc

format short e
addpath(genpath(pwd))

E = 1; C = 1; R = 1;

% Test
PDE.term1.Opt = 'Grad';
PDE.term1.FunCoef = @(x)(1-x.^2);
PDE.term1.Coef = -E;

PDE.term2.Opt = 'Diff';
PDE.term2.FunCoef = @(x)(1-x.^2);
PDE.term2.Coef = C;

PDE.term3.Opt = 'Grad';
PDE.term3.FunCoef = @(x)(x.*(1-x.^2));
PDE.term3.Coef = -R;



Lev = 4;
Deg = 2;
num_plot = 3;

LInt = -1;
LEnd = 1;
% Lmax = Lend-Lstart;

%% Matrix
% Term 1
Mat_Term1 = MatrixGrad(Lev,Deg,LInt,LEnd,PDE.term1.FunCoef);
Mat_Term2 = MatrixDiff(Lev,Deg,LInt,LEnd,PDE.term2.FunCoef,0);
Mat_Term3 = MatrixGrad(Lev,Deg,LInt,LEnd,PDE.term3.FunCoef);
Mat = PDE.term1.Coef*Mat_Term1 + ...
    PDE.term2.Coef*Mat_Term2 + ...
    PDE.term3.Coef*Mat_Term3 ;
%% RHS

%% B.C


%% Solve

%--Quadrature
quad_num=10;
%---------------

% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre(1,Deg);

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);

%---------------------------
% Jacobi of variable x and v
% Define Matrices
%---------------------------


n=2^(Lev);h=Lmax/n;
Jacobi=h;
dof_1D=Deg*n;
A12 = sparse(dof_1D,dof_1D);
f0 = sparse(dof_1D,1);

b = sparse(dof_1D,1);bb = sparse(dof_1D,1);
fexact = sparse(dof_1D,1);
qexact = sparse(dof_1D,1);

CFL = 0.001;
% dt = CFL*h^((Deg-1)/3)/2;
EndTime = 3;
dt = CFL*h^((Deg)/3);
maxT = ceil(EndTime/dt)

% Assume
% [ I  A12]
% [A21  0 ]
% as the matrix
% e.g. we assume neuman boundary:: q=0 on boundary
%
B12 = sparse(dof_1D,dof_1D);
A21 = sparse(dof_1D,dof_1D);
A12 = sparse(dof_1D,dof_1D);
% generate 1D matrix for DG
for L=0:n-1
    
    %---------------------------------------------
    % (funcCoef*q,d/dx p)
    %---------------------------------------------
    x0 = Lstart+L*h;
    x1 = x0+h;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;
    xmid= (x0+x1)/2;
    
    val=1/h*[Dp_val'*(quad_w.*funcCoef(xi).*p_val)];
    
    c = Deg*L+1:Deg*(L+1);
    
    B12 = B12 + sparse(c'*ones(1,Deg),ones(Deg,1)*c,val,dof_1D,dof_1D);
    
    % diffusion term
    val=1/h*[Dp_val'*(quad_w.*funcCoef(xi ).*p_val)]+...
              p_val'*(quad_w.*funcCoef2(xi).*p_val)/2;    
    A12 = A12 + sparse(c'*ones(1,Deg),ones(Deg,1)*c,val,dof_1D,dof_1D);
    
    val=1/h*[Dp_val'*(quad_w.*p_val)];
    A21 = A21 + sparse(c'*ones(1,Deg),ones(Deg,1)*c,val,dof_1D,dof_1D);
    
    
    val = sqrt(h)/2*[p_val'*(quad_w.*exactf(xi,0))];
    fexact(c)=fexact(c)+val;
    
    
    %----------------------------------------------
    % -<funcCoef*{q},p>
    %----------------------------------------------
    if FluxType == 'CF'
%         val=[ p_1'*funcCoef(x0)*p_2   p_1'*funcCoef(x0)*p_1,...
%             -p_2'*funcCoef(x1)*p_2  -p_2'*funcCoef(x1)*p_1]/2/h;
%         B12=B12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
%             val(:,Deg+1:2*Deg)+val(:,2*Deg+1:3*Deg),...
%             dof_1D,dof_1D);
        
        val=[p_1'*funcCoef(x0)*p_2   -p_2'*funcCoef(x1)*p_2]/h;
        B12=B12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
            val(:,Deg+1:2*Deg),...
            dof_1D,dof_1D);
        
        val2=[ p_1'*funcCoef(x0)*p_2   p_1'*funcCoef(x0)*p_1,...
            -p_2'*funcCoef(x1)*p_2  -p_2'*funcCoef(x1)*p_1]/2/h;
        A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
            val2(:,Deg+1:2*Deg)+val2(:,2*Deg+1:3*Deg),...
            dof_1D,dof_1D);
        
        val3 = [ p_1'*p_2   p_1'*p_1,...
            -p_2'*p_2  -p_2'*p_1]/2/h;
        A21=A21+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
            val3(:,Deg+1:2*Deg)+val3(:,2*Deg+1:3*Deg),...
            dof_1D,dof_1D);
        
        
        if L>0
%             B12 = B12+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val(:,1:Deg),dof_1D,dof_1D);
            B12=B12+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val(:,1:Deg),dof_1D,dof_1D);
                        
            A12 = A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val2(:,1:Deg),dof_1D,dof_1D);
            A21 = A21+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val3(:,1:Deg),dof_1D,dof_1D);
        elseif L == 0
            % new implementation
            A21 = A21 +sparse(c'*ones(1,Deg),ones(Deg,1)*c,-p_1'*p_1/h/2,dof_1D,dof_1D);
        end
        if L<n-1
%             B12 = B12+sparse(c'*ones(1,Deg),ones(Deg,1)*c+Deg,val(:,3*Deg+1:4*Deg),dof_1D,dof_1D);
            
            A12 = A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c+Deg,val2(:,3*Deg+1:4*Deg),dof_1D,dof_1D);
            A21 = A21 +sparse(c'*ones(1,Deg),ones(Deg,1)*c,-p_1'*p_1/h/2,dof_1D,dof_1D);
        elseif L == n-1
            
            
            A21=A21+sparse(c'*ones(1,Deg),ones(Deg,1)*c,p_2'*p_2/h/2,dof_1D,dof_1D);
        end
    end
    
    % up-winding flux
    if FluxType == 'UF'

        val=[p_1'*funcCoef(x0)*p_2   -p_2'*funcCoef(x1)*p_2]/h;
        B12=B12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
            val(:,Deg+1:2*Deg),...
            dof_1D,dof_1D);
        
        val2=[p_1'*funcCoef(x0)*p_1   -p_2'*funcCoef(x1)*p_1]/h;
        A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
            val2(:,1:Deg),...
            dof_1D,dof_1D);

        val3 = [ p_1'*p_2 -p_2'*p_2]/h;
        A21=A21+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
            val3(:,Deg+1:2*Deg),...
            dof_1D,dof_1D);
        
        if L>0
            B12=B12+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val(:,1:Deg),dof_1D,dof_1D);
            %A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val(:,1:Deg),dof_1D,dof_1D);
            A21=A21+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val3(:,1:Deg),dof_1D,dof_1D);
        end
        if L<n-1
            A12 = A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c+Deg,val2(:,1+Deg:2*Deg),dof_1D,dof_1D);
        end
    end
    
    
    
    val = sqrt(h)/2*[p_val'*(quad_w.*ff0(xi))]; %exactf(xi,0))];
    f0(c) = val;
    
    
end



[quad_x,quad_w]=lgwt(num_plot,-1,1);
% quad_x = [-1,1]';

p_val = legendre(quad_x,Deg);
for L=0:n-1
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    
    Iu = [Deg*L+1:Deg*(L+1)];
    
    Iv = [num_plot*L+1:num_plot*(L+1)];
    %     xi=h*(quad_x/2+1/2+L);
    
    x0 = Lstart+L*h;
    x1 = x0+h;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;[L*h,L*h+h];%
    
    
    Meval(Iv,Iu)=sqrt(1/h)*p_val;
    x_node(Iv,1)=xi;
    
end

% checked of projection
plot(x_node,Meval*f0,'r-o',x_node,Meval*(A12*f0),'b-o',x_node,exactf(x_node,0),'b--','LineWidth',2);

% return
Mat = E*B12+C*A21*A12;

% max(abs(Meval*f0))
% val = Meval*A12*f0-exactq(x_node,0);
% [norm(val) max(abs(val))]

b = b+bb;


[quad_x,quad_w]=lgwt(num_plot,-1,1);
total_particle = 0;
ffval = Meval*f0;

for i = 1:num_plot
    total_particle =  total_particle+quad_w(i)*h/2*sum(ffval(i:num_plot:end));
end
total_particle

tp(1) = total_particle;

figure
for t = 1:maxT
    %     t
    time = t*dt;
    %     tmp = A12'*A12*f0;%dt*A12'*A12*f0+dt*b*exp(time);
    %     tmp2 = A12*A12*f0;
    
    % %     fval = f0+dt*Mat*f0+dt*b;%*exp(time);
    
    f1 = f0 + dt*( Mat*f0+b*exp(time-dt) );
    f2 = 3/4*f0+1/4*f1+1/4*dt*(Mat*f1+b*exp(time));
    fval = 1/3*f0+2/3*f2+2/3*dt*(Mat*f2+b*exp(time-dt/2));
    
    %     f1 = f0 + dt*( Mat*f0+b );
    %     f2 = 3/4*f0+1/4*f1+1/4*dt*(Mat*f1+b );
    %     fval = 1/3*f0+2/3*f2+2/3*dt*(Mat*f2+b );
    
    f0 = fval;
    
    plot(x_node,Meval*f0,'r-o',x_node,Meval*(A12*f0),'b-<',x_node,exactf(x_node,time),'r--');
%     hold on
%     plot(x_node,exactf(x_node,time),'r--')
    title(['time at ',num2str(time)])
    pause (0.1)
    
    total_particle = 0;
    ffval = Meval*f0;
    for i = 1:num_plot
        
        total_particle =  total_particle+...
            quad_w(i)*h/2*sum(ffval(i:num_plot:end));
    end
    tp(t+1) = total_particle;
    
    if abs(time-0.5)<=dt || abs(time-1)<=dt || abs(time-1.5)<=dt || abs(time-2)<=dt || abs(time-2.5)<=dt ||abs(time-3)<=dt
        save(['hyper_',FluxType,'_Deg',num2str(Deg),'_Lev',num2str(Lev),'_End',num2str(time),'.mat'])
    end
end
% figure;plot(x,f_loc'*f0,'r-o');hold on;
% plot(x,exactf(x,time),'b--')
hold on
plot(x_node,exactf(x_node,time),'r-o')
%  max(abs(Meval*f0))
val = Meval*f0-exactf(x_node,time);
%  [norm(val) max(abs(val))]

fL2 = 0; fLinf = max(abs(val));
%  total_particle = 0;
ffval = Meval*f0;
for i = 1:num_plot
    fL2 = fL2 + quad_w(i)*h/2*sum(val(i:num_plot:end).^2);
    %       total_particle =  total_particle+quad_w(i)*h/2*sum(ffval(i:num_plot:end));
end
[sqrt(fL2) fLinf]
%  total_particle


%  err = f0-fexact;
%  full([norm(err) max(abs(err))])
%
%  err = A12*f0-qexact;
%  full([norm(err) max(abs(err))])

figure;
plot(x_node,Meval*f0,'r-o',x_node,exactf(x_node,time),'r--','LineWidth',2);


figure;
plot(tp,'r-o')
