%%% UserFile
% This is the main interface file for HyPaD. Please make sure that you have
% installed and started Gurobi before calling HyPaD. If you want that the
% initial box B = [z,Z] is computed by interval arithmetic, you also need
% to install and start INTLAB before calling HyPaD.
% Call initSession.m prior to first use

%% Some clean-up first
clear;
close all;
clc;

%% Please enter your parameters below
% Your Problem
problem = 'T4';
param = [4;4];

% Provide initial bounds yourself or leave empty to auto-compute (INTLAB is
% required for auto-compute)
L = [];
U = [];

% Set quality epsilon and offset
EPSILON = 0.1;
OFFSET = EPSILON*1e-3;

% Specify [mode, options, guess] for SNIA in HyPaD
% Mode [1 == fixed boxes, 2 == dynamic boxes, 3 == list]
% Options [number of boxes for fixed boxes, 
%          point for dynamic boxes: 1 == lb, 2 == ub, 3== mid,  
%          no options for list]
% Guess [0 == none, 1 == based on continuous relaxation]
HYPAD_PARAM = [1,4,0];

% Should the result be plotted (m = 2 and m = 3 only) [1 == yes, 0 == no]
plot_result = 1;

%% Call solver
[L,U,N,ids,it,exitflag,time] = callSolver(problem,param,L,U,HYPAD_PARAM,EPSILON,OFFSET,plot_result)