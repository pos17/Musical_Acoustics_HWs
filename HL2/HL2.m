clear; close all; clc;
% Set of parameters

V0 = 0.1; % m3
l = 10e-2; % m
S = 100e-6; % m2
c = 343; % m/s
rho = 1.2; % kg/m3

Fs = 48000;
dur = 3;
N = dur*Fs+1;

%% EX1 - electrical

M =  rho*l/S;
C = V0/(rho*c^2);
R = rho*c/S;

set_param('Ex1', 'PreLoadFcn', num2str(Fs))
set_param('Ex1/M1', 'l', num2str(M));
set_param('Ex1/C1', 'c', num2str(C));
set_param('Ex1/R1', 'R', num2str(R));


open_system("HL2\Ex1.slx", 'loadonly');
out = sim("HL2\Ex1.slx", dur);

input = out.force.Data;
output = out.velocity.Data;

f = linspace(0, Fs, N);
H = fft(output)./fft(input);
f0 = f(find(db(abs(H))==max(db(abs(H))),1));

close all
figure('Renderer', 'painters', 'Position', [100 100 800 400])
plot(f, db(abs(H)), LineWidth=1.2);
xlim([0, 2*f0]);% ylim([-100, 0]);
xlabel('Freq [Hz]'); ylabel("|H| [dB]");
title("Frequency response of the resonator")
hold on
xline(f0, '--')
text(f0*(1.01), -70, 'f_0')
grid minor
% delete(".\plots\Ex1_A_FRF.png");
% saveas(gcf, ".\plots\Ex1_A_FRF.png");

%% EX1 B

r = sqrt(S/pi);
deltaL = 1.6*r; 
f0_an = c/(2*pi)*sqrt(S/(l*V0));

