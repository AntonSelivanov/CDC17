% This MATLAB program checks the feasibility of LMIs from Proposition 1 of the paper 
% A. Selivanov and E. Fridman, "Delayed boundary control of a heat equation under discrete-time point measurements," 
% in 56th IEEE Conference on Decision and Control, 2017, pp. 1248–1253.

a=10;                       % reaction coefficient in (1)
dL=1; dR=0;                 % boundary conditions: Dirichlet on the left, Neumann on the right 
r=.05;                      % input delay 
h=.01;                      % sampling period 
etaM=.01;                   % maximum network delay
OmegaM=1/10;                % maximum subdomain size for N=10 sensors
L=-10;                      % injection gain in (4)
alpha0=.48;                 % decay rate of the observation error in (5)
alpha1=1;                   % tuning parameter for Halanay inequality 

if LMI_CDC17_prop1(a,dL,dR,r,h,etaM,OmegaM,L,alpha0,alpha1)
    disp('Proposition 1: Feasible'); 
else
    disp('Proposition 1: Not feasible'); 
end
