clc; clear all; close all;
%% When use option 1, comment all contents for option 2. When using option 2, 
%% comment all contents for option 1
%% Option 1: define all quantities with unit
E = 1.8e6;
Poisson = 0.5;
h = 1.6e-3;
rho = 1180;
g = 10;
Lgb = (h^2 * E/(8*rho*g))^(1/3);
ls = 30;

ks = E * pi *h^2;
kb = E * pi * h^4/4;
normks = 200;
kappaT = 10 * Lgb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Option 2: Use nondimensionless group
% User defined variables
% Poisson = 0.5;
% normks = 1000;
% kappaT = 1;
% ls = 15;
% 
% Lgb  = 1;
% ks = normks;
% kb = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xTop, yTop, zTop, alpha] = getLocalOpt(ls, 0, Lgb, ks, kb);

fprintf('hanglength %f, norm ks %f, target curvature %f, Poisson ratio %f with intial guess (%f, %f, %f)\n', ...
          ls, normks, kappaT, Poisson, xTop, yTop, zTop);
