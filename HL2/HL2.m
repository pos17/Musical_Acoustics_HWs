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

r = sqrt(S/pi);
deltaL = 0.6*r+(8/(3*pi)*r); 
l1 = l+deltaL;  
% l1=l;

f = linspace(0, Fs, N);

%% EX1 - electrical

M =  rho*l1/S;
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
ylim([-100, 0])
% delete(".\plots\Ex1_A_FRF.png");
% saveas(gcf, ".\plots\Ex1_A_FRF.png");

%% EX1 B


f0_an = c/(2*pi)*sqrt(S/(l1*V0));
% Z0 = rho*c/S;
% w = 2*pi*f;
% nu = c*(1-(1.65e-3)./(r*sqrt(f)));
% alpha = 3e-5*sqrt(f)/r;
% % num = tanh(alpha*l)+1i*tan(w*l./nu);
% % den = 1 + 1i*tanh(alpha*l).*tan(w*l./nu);
% 
% k = w/c;
% Z_v = rho*c^2./(1i*w*V0);
% 
% num = Z_v.*cos(k*l) + 1i*Z0*sin(k*l);
% den = 1i*Z_v.*sin(k*l) + Z0*cos(k*l);
% 
% Zin = Z0*num./den;

