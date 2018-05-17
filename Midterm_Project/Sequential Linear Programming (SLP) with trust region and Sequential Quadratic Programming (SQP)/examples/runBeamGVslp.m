%% runBeamGVslp
% Script to Run Gary Vanderplaats cantilever Beam with slp_trust. 
% N beam segments and 2N design variables, taken from 
% Vanderplaats (1984) Example 5-1, pp. 147-150.
%
%  Trust Region (TR) Strategy follows
%  Nocedal, J., and Wright, S.J. Numerical Optimization. 
%  New York: Springer, 2006.
%  Algorithm  4.1 for TR and
%  Algorithm 15.1 for filter; and 
%  Equation (15.4) for the simple L1 merit function
%
% TR Approximate Model (TRAM) for adaptive move limits roughly follows 
% Wujek, B. A., and Renaud, J. E. "New Adaptive Move-Limit Management
% Strategy for Approximate Optimization, Part I," AIAA J. Vol. 36, No. 10,
% 1998, pp. 1911-1921 
% except that quadratic approximations are used for objective and
% constraints and a linear approximation for the merit penalty function
% (instead of TANA for the merit function)
% during the "line search" that determines the trust region radius.
% Also, a custom TR algorithm is used for multi-penalty L1 merit function. 

%% Initialize variables
clear; clc
N = 10;%[5 10 20 40 50 100 200]; % number of beam segments
x0 = [5*ones(1,N), 40*ones(1,N)];
xlb = [ ones(1,N),  5*ones(1,N)];

%% SLP Trust
%  TrustRegion='filter' (default)
options.Display = 'Iter';
options.TolX = .5;
options.TolFun = 0.5;
options.TolCon = 1e-4;
options.MoveLimit = 0.2;
options.MaxIter = 40;
% default options for simple merit function TR algorithm
[dv1,f1,sta1,out1,lambda1] = slp_trust(@fbeamGV,x0,options,xlb,[],@gbeamGV);

%% other SLP TR methods - simple merit descent function
%  Exact L1 single penalty, mu*max(0,g) (default)
options.TrustRegion = 'simple';
options.Contract = 0.2;
slp_trust(@fbeamGV,x0,options,xlb,[],@gbeamGV);

%% TrustRegion='merit' - multiple-penalty parameter L1 merit function
options.TrustRegion = 'merit';
slp_trust(@fbeamGV,x0,options,xlb,[],@gbeamGV);

%% TrustRegion='TRAM' - TR adaptive move limit bounds (with merit)
options.TrustRegion = 'TRAM';
slp_trust(@fbeamGV,x0,options,xlb,[],@gbeamGV);

%% TrustRegion='TRAM-filter' combination
options.TrustRegion = 'TRAMfilter';
slp_trust(@fbeamGV,x0,options,xlb,[],@gbeamGV);