%%% UserFile
% This is the main interface file for HyPaD. Please make sure that you have
% installed and started Gurobi before calling HyPaD. If you want that the
% initial box B = [z,Z] is computed by interval arithmetic, you also need
% to install and start INTLAB before calling HyPaD.

%% Some clean-up first
clear;
close all;
clc;

%% Please enter your parameters below
% Your Problem
problem = 'H1';
param = [2;4];

% Set quality epsilon and offset
EPSILON = 0.1;
OFFSET = EPSILON*1e-3;

% Provide initial bounds yourself or leave empty to auto-compute (INTLAB is
% required for auto-compute)
z = [];
Z = [];

% Should the result be plotted (m = 2 and m = 3 only) [1 == yes, 0 == no]
plot_result = 1;

%% Call solver
[L,U,N,ids,it,flag,time,z,Z] = callSolver(problem,param,z,Z,EPSILON,OFFSET,plot_result);