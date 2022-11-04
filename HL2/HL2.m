clear; close all; clc;
% Set of parameters

V0 = 0.1; % m3
l = 10e-2; % m
S = 100; % m2
c = 343; % m/s
rho = 1.2; % kg/m3

Fs = 48000;
dur = 3;
N = dur*Fs+1;

%% EX1 - electrical

M =  rho*l/S;
C = V0/(rho*c^2);

set_param('Ex1/M1', 'l', num2str(M));
set_param('Ex1/C1', 'c', num2str(C));


open_system("HL2\Ex1.slx", 'loadonly');
out = sim("HL2\Ex1.slx", dur);
%%
input = out.force.Data;
output = out.velocity.Data;

f = linspace(0, Fs, N);
H = fft(output)./fft(input);

close all
figure;
plot(f, db(abs(H)), LineWidth=1.2);
xlim([0, Fs/2]); ylim([-100, 0]);
xlabel('Freq [Hz]'); ylabel("|H| [dB]");
title("")
grid minor



