%This is for plotting the error between the error between exp(A*dt)*v and
%IIF exp

lev = 7; %level
T = 5*1/2^(lev);
deg = 4; %degree
interval = 5;
num = 60;%number of values of M tested
Mmax = interval*num; %maximum M

%dt = 2;
dt = 1/2^(lev); %time step
n = 5; %Number of time steps
%n = ceil(T/dt);

%combined electric field and collision
%asgard(fokkerplanck1_4p4,'lev',lev,'deg',deg,'implicit',true,'num_steps',1)

%just electric field
%asgard(fokkerplanck1_4p1a,'lev',lev,'deg',deg,'implicit',true,'num_steps',1)

%just collision
%asgard(fokkerplanck1_4p2,'lev',lev,'deg',deg,'implicit',true,'num_steps',1)

asgard(mirror_velocity2,'timestep_method','BE', 'dt', 1e-5, 'num_steps',1, 'grid_type', 'SG', 'deg', 4, 'lev', 5)

f0 = load('initial_vec000.mat','f0');
f0 = f0.f0; %f0 is initial vector

load('pde.mat','pde')
load('opts.mat','opts')
load('hash_table.mat','hash_table')
load('Meval.mat','Meval')
load('nodes.mat','nodes')
load('matrix_iter000.mat','A');

AN = expm(A*T)*f0;
an = wavelet_to_realspace(pde,opts,Meval,AN,hash_table);
M = linspace(interval,Mmax,num);
E = zeros(1,num);

[N,~] = size(A);
F = zeros(N,num);

for j=1:num
    m = interval*j;
    %IIF
    f = f0;
    for i=0:n-1
        [V,H] = myarnoldi(A,f,m);
        gamma = norm(f);
    
        kryExp = expm(H*dt);
        kryExp = kryExp(:,1); %first column in matrix
        kryExp = gamma*V*kryExp;

        f = kryExp;
    end
    f1 = wavelet_to_realspace(pde,opts,Meval,f,hash_table); %Krylov 2nd Order
    F(:,j) = f1;
    E(j) = norm(f1-an)/sqrt(N);
end

semilogy(M,E);