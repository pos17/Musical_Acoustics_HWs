%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HL2
% Exercise 1
% Helmoholtz resonators tree
% 
% This script calls the simulink file that implements the simple electrical
% equivalent of a Helmholtz resonator.
%
% OPEN SIMULINK MODEL Ex1.slx AND THEN RUN THIS SCRIPT
% 
% 
% Musical Acoustic Course
% Donà Stefano
% Ostan Paolo
% 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
% Set of parameters
V0 = 0.1; % m3
l = 10e-2; % m
S = 100; % m2
c = 343; % m/s
rho = 1.2; % kg/m3

% virtual elongation of the neck
r = sqrt(S/pi);
deltaL = 0.6*r+(8/(3*pi)*r); 
l1 = l+deltaL;  

% sampling frequency and duration of signal 
Fs = 48000;
dur = 3;
N = dur*Fs+1;

% frequency axis
f = linspace(0, Fs, N);
w = 2*pi*f;

%% EX1 A - electrical simulation

% electrical equivalents
M =  rho*l1/S;
C = V0/(rho*c^2);
R = rho*c/S;

% setting parameters in the model
set_param('Ex1', 'PreLoadFcn', num2str(Fs))
set_param('Ex1/L1', 'l', num2str(M));
set_param('Ex1/C1', 'c', num2str(C));
set_param('Ex1/R1', 'R', num2str(R));

% simulation
open_system("HL2\Ex1.slx", 'loadonly');
out = sim("HL2\Ex1.slx", dur);

% extracting data from the simulation results
input = out.force.Data;
output = out.velocity.Data;

% mobility computation
H = fft(output)./fft(input);
f0 = f(find(db(abs(H))==max(db(abs(H))),1));

% numerical admittance computation
Z0 = R + 1i*(w*M-1./(w*C));
Y0 = 1./Z0;

%% PLOTTING
close all
figure('Renderer', 'painters', 'Position', [100 100 800 400])
plot(f, db(abs(H)), LineWidth=1.7);
xlim([0, 2*f0]);
xlabel('Freq [Hz]'); ylabel("|H| [dB]");
title("Frequency response of the resonator")
hold on
xline(f0, '--', LineWidth=1.2)
text(f0*(1.01), -70, 'f_0', FontSize=12)
grid minor
ylim([-100, 0])
hold on
plot(f, db(abs(Y0)), 'r--', LineWidth=1.8);
legend('experimental', 'numerical')
% delete(".\plots\Ex1_A_FRF.png");
% saveas(gcf, ".\plots\Ex1_A_FRF.png");

%% EX1 B - numerical computation of resonant frequency

f0_an = c/(2*pi)*sqrt(S/(l1*V0));

