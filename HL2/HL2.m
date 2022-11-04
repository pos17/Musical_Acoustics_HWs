clear; close all; clc;
% Set of parameters

V0 = 0.1; % m3
l = 10e-2; % m
S = 100; % m2
c = 343; % m/s
rho = 1.2; % kg/m3

Fs = 48000;
dur = 3;
N = dur*Fs;

%% EX1 - electrical

M =  rho*l/S;
C = V0/(rho*c^2);