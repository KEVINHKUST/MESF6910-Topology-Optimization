%% runSpringExampleTrust.m
% Example problem taken from Vanderplaats textbook, example 3-1. 
% Unconstrained potential energy minimization of two springs.
% Complex-step gradient.

%% Initialize tolerances
clear; clc;
options.ComplexStep = 'on';
options.Display = 'iter';
options.MaxIter = 50;
options.TypicalX = [5;5];
options.TolFun  = 0.001; % 0.0001 for SLP slow termination criterion
options.TolX    = 0.1;

%% SQP
[x,out]=sqp(      @fVanderplaatsSpringEx3d1,[5 5],options)

%% SLP gradually reduced move limits (no Trust Region Strategy)
options.TrustRegion = 'off';
slp_trust(@fVanderplaatsSpringEx3d1,[5 5],options)

%% SLP Trust Region
options.TrustRegion = 'TRAM'; % Trust Region Adaptive Move Limits
options.MoveLimit   = 0.5;
options.OptimalityTolerance = 0.1;
[X,PE] =slp_trust(@fVanderplaatsSpringEx3d1,[5 5],options)