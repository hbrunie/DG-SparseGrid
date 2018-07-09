%===================================================
% This code solves Maxwell's equation by
% DG-SG grids methods
%===================================================
%clc
clear
close all

format short e

addpath(genpath(pwd))

%% Step 1. Setting Parameters
Lev = 3;
Deg = 2;


Lmax = 1;

pde = Maxwell1;
MaxT = 100;
dt = 1/10000;


%*************************************************
%% Step 1.1. Set Up Matrices for Multi-wavelet
% Input:  Deg and Lev
% Output: Convert Matrix FMWT_COMP
%*************************************************
FMWT_COMP_x = OperatorTwoScale(Deg,2^Lev);

%*************************************************
%% Step 1.2. Set Up Initial Conditions
% F_1D: rhs coefficients
% E_1D: Initial Condition for E
% B_1D: Initial condition for B
% Input:  Deg, Lev,Lmax, pde,FMWT_COMP_x
% Output: F_1D,E_1D,B_1D
%*************************************************
[F_1D,E_1D,B_1D] = Intial_Con(Lev,Deg,Lmax,pde,FMWT_COMP_x);

%*************************************************
%% Step 2. Hash Table and 1D Connectivity
%*************************************************
[Hash,InvHash]=HashTable(Lev);
%1D connectivity
Con1D=Connect1D(Lev);

%*************************************************
%% Step 3. Coefficient Matrix for Time-independent Matrix
%*************************************************
GradX = Matrix_TI(Lev,Deg,Lmax,FMWT_COMP_x);

%*************************************************
%% Step 4.3D Maxwell Solver
%*************************************************
% global rhs, E0, and B0 vectors
% b_s is the RHS vector
% E_s and B_s are used for error estimate

%% Maxwell Solver
[Bh,E1h,E2h] = MaxwellSolver(Lev,Deg,Hash,InvHash,Con1D,GradX,pde.w,dt,MaxT,...
    F_1D,E_1D.E1*cos(0),E_1D.E2*cos(0),B_1D.B*0);
sol_n=[Bh;E1h;E2h];

%% Error Estimate
time=dt*MaxT;
u_s=[B_1D.B*sin(pde.w*time);E_1D.E1*cos(pde.w*time);E_1D.E2*cos(pde.w*time)];

full([Deg Lev max(abs(sol_n-u_s)) norm(sol_n-u_s)])

