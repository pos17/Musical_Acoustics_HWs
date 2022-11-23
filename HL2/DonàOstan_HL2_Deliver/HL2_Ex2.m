%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HL2
% Exercise 2
% Helmoholtz resonators tree
% All the parameters assigned are imported from the exercise 1
% This script call the simulink file that implements the Helmotz trees
% discussed in the report, the trees are not in the same order but the
% topology and parameters are specified by the comments on the code and 
% on the graphs.
% Before the run of the matlab code please open on simulink the files: 
%   * Ex2_A.slx
%   * Ex2_B.slx
%   * Ex2_C.slx
%   * Ex2_D.slx
%   * Ex2_E.slx
%
% 
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
l1 = 10e-2; % m
S = 100; % m2
c = 343; % m/s
rho = 1.2; % kg/m3

Fs = 48000;
dur = 3;
N = dur*Fs+1;
a = sqrt(S/pi);

l = l1 + 8/3/pi*a + 0.6*a;
%l=l1
%% Electrical components 

M =  rho*l/S;
C = V0/(rho*c^2);
R = rho*c/S;

f = linspace(0, Fs, N);
w = 2*pi*f;

Z = R+1i*w*M+(1i*w*C).^(-1);

%% EX2 A Full 3x2 Tree
set_param('Ex2_A', 'PreLoadFcn', num2str(Fs));


for i=1:7
    subsys_name =strcat('Ex2_A/Subsystem',num2str(i),'/');
    set_param(strcat(subsys_name,'L1'), 'l', num2str(M));
    set_param(strcat(subsys_name,'C1'), 'c', num2str(C));
    set_param(strcat(subsys_name,'R1'), 'R', num2str(R));
end

open_system("Ex2_A.slx", 'loadonly');
out = sim("Ex2_A.slx", dur);
disp("done")

vel0 = out.vel0.Data;
vel1 = out.vel1.Data;
vel2 = out.vel2.Data;
press = out.pressure.Data;

H0 = db(abs(fft(vel0)./fft(press)));
H1 = db(abs(fft(vel1)./fft(press)));
H2 = db(abs(fft(vel2)./fft(press)));

[val,f0] = findpeaks(H0);
f0 = 48000/144001 * f0;
%close all
figure('Renderer', 'painters', 'Position', [100 100 1000 500])
plot(f, H0, LineWidth=1.2);
xlim([0, 2500]); ylim([-100, 0])
hold on

plot(f, H1, LineWidth=1.2);
hold on
plot(f, H2, LineWidth=1.2);
grid minor
hold on
xline(f0(1), 'k--')
text(f0(1)*(1.01), -90, strcat('f_0=',num2str(f0(1))))
xline(f0(2), 'k--')
text(f0(2)*(1.01), -90, strcat('f_1=',num2str(f0(2))))
xline(f0(3), 'k--')
text(f0(3)*(1.01), -90, strcat('f_2=',num2str(f0(3))))
legend('H0', 'H1', 'H2')
xlabel('Freq [Hz]'); ylabel("|H| [dB]");
title("Full 3x2 Tree, Frequency Response functions of the system")

%delete("./plots/Ex2_A_FRF.png");
%saveas(gcf, "./plots/Ex2_A_FRF.png");

%% EX2 B 3x3 Tree 

set_param('Ex2_B', 'PreLoadFcn', num2str(Fs))
for i=1:7
    subsys_name =strcat('Ex2_B/Subsystem',num2str(i-1),'/');
    set_param(strcat(subsys_name,'L1'), 'l', num2str(M));
    set_param(strcat(subsys_name,'C1'), 'c', num2str(C));
    set_param(strcat(subsys_name,'R1'), 'R', num2str(R));
end

open_system("Ex2_B.slx", 'loadonly');
out2_B = sim("Ex2_B.slx", dur);

press_B = out2_B.pressure.Data;
vel0_B = out2_B.vel0.Data;
vel1_B = out2_B.vel1.Data;
vel2_B = out2_B.vel2.Data;
vel3_B = out2_B.vel3.Data;

f = linspace(0, Fs, N);
H0 = db(abs(fft(vel0_B)./fft(press_B)));
H1 = db(abs(fft(vel1_B)./fft(press_B)));
H2 = db(abs(fft(vel2_B)./fft(press_B)));
H3 = db(abs(fft(vel3_B)./fft(press_B)));

%close all;
[val,f0] = findpeaks(H0);
f0 = 48000/144001 * f0;
figure('Renderer', 'painters', 'Position', [100 100 1000 400])
plot(f, H0, LineWidth=1.2);
xlim([0, 2500]); ylim([-110, 0]);
hold on
plot(f, H1, LineWidth=1.2);
hold on
plot(f, H2, LineWidth=1.2);
hold on
plot(f, H3, LineWidth=1.2);
xlabel('Freq [Hz]'); ylabel("|H| [dB]");
xline(f0(1), 'k--')
text(f0(1)*(1.01), -90, strcat('f_0=',num2str(f0(1))))
xline(f0(2), 'k--')
text(f0(2)*(1.01), -90, strcat('f_1=',num2str(f0(2))))
xline(f0(3), 'k--')
text(f0(3)*(1.01), -90, strcat('f_2=',num2str(f0(3))))
xline(f0(4), 'k--')
text(f0(4)*(1.01), -90, strcat('f_3=',num2str(f0(4))))
%text(f0*(1.01), -70, 'f_0')
grid minor
legend('H0', 'H1', 'H2', 'H3')
title("3x3 Tree, Frequency Response functions of the system")
%delete(".\plots\Ex2_B_FRF.png");
%saveas(gcf, ".\plots\Ex2_B_FRF.png");


%% EX2 C 2x3 Tree 

set_param('Ex2_C', 'PreLoadFcn', num2str(Fs))
for i=1:4
    subsys_name =strcat('Ex2_C/Subsystem',num2str(i-1),'/');
    set_param(strcat(subsys_name,'L1'), 'l', num2str(M));
    set_param(strcat(subsys_name,'C1'), 'c', num2str(C));
    set_param(strcat(subsys_name,'R1'), 'R', num2str(R));
end

open_system("Ex2_C.slx", 'loadonly');
out2_C = sim("Ex2_C.slx", dur);

press_C = out2_C.pressure.Data;
vel0_C = out2_C.vel0.Data;
vel1_C = out2_C.vel1.Data;

f = linspace(0, Fs, N);
H0 = db(abs(fft(vel0_C)./fft(press_C)));
H1 = db(abs(fft(vel1_C)./fft(press_C)));


%close all;

[val,f0] = findpeaks(H0);
f0 = 48000/144001 * f0;

figure('Renderer', 'painters', 'Position', [100 100 1000 400])
plot(f, H0, LineWidth=1.2);
xlim([0, 2500]); ylim([-110, 0]);
hold on
plot(f, H1, LineWidth=1.2);
xlabel('Freq [Hz]'); ylabel("|H| [dB]");
xline(f0(1), 'k--')
text(f0(1)*(1.01), -90, strcat('f_0=',num2str(f0(1))))
xline(f0(2), 'k--')
text(f0(2)*(1.01), -90, strcat('f_1=',num2str(f0(2))))
grid minor
legend('H0', 'H1')
title("2x3 Tree, Frequency Response functions of the system")
%delete(".\plots\Ex2_C_FRF.png");
%saveas(gcf, ".\plots\Ex2_C_FRF.png");

%% EX2 D Full 2x2 Tree 


% Numerical Solution for tree of height 2 and two siblings every parent

Z_sibs = ((((1i*2*pi*f).^2)*M*C)+(1i*2*pi*f*R*C)+1)./(1i*2*pi*f*C);


Z_2 = (1i*2*pi*f*M) + R +(1./((1i*2*pi*f*C) + (2./(Z_sibs))));
H2 =  db(1./abs(Z_2));


set_param('Ex2_D', 'PreLoadFcn', num2str(Fs))
for i=1:3
    subsys_name =strcat('Ex2_D/Subsystem',num2str(i-1),'/');
    set_param(strcat(subsys_name,'L1'), 'l', num2str(M));
    set_param(strcat(subsys_name,'C1'), 'c', num2str(C));
    set_param(strcat(subsys_name,'R1'), 'R', num2str(R));
end

open_system("Ex2_D.slx", 'loadonly');
out2_D = sim("Ex2_D.slx", dur);

press_D = out2_D.pressure.Data;
vel0_D = out2_D.vel0.Data;
vel1_D = out2_D.vel1.Data;

f = linspace(0, Fs, N);
H0 = db(abs(fft(vel0_D)./fft(press_D)));
H1 = db(abs(fft(vel1_D)./fft(press_D)));


%close all;

[val,f0] = findpeaks(H0);
f0 = 48000/144001 * f0;

figure('Renderer', 'painters', 'Position', [100 100 1000 400])
plot(f, H0, LineWidth=1.2);
xlim([0, 2500]); ylim([-110, 0]);
hold on
plot(f, H1, LineWidth=1.2);
% xlim([0, Fs/2]);

hold on
plot(f, H2, LineWidth=1.2,LineStyle="-.");
% xlim([0, Fs/2]);
hold  on
xlabel('Freq [Hz]'); ylabel("|H| [dB]");
%xline(f0, '--')
xline(f0(1), 'k--')
text(f0(1)*(1.01), -90, strcat('f_0=',num2str(f0(1))))

xline(f0(2), 'k--')
text(f0(2)*(1.01), -90, strcat('f_1=',num2str(f0(2))))
grid minor
legend('H0', 'H1', 'H0 numerical')
title("Full 2x2 Tree, Frequency Response functions of the system")
%delete(".\plots\Ex2_D_FRF.png");
%saveas(gcf, ".\plots\Ex2_D_FRF.png");

%% EX2 E 4x2 Tree 


set_param('Ex2_E', 'PreLoadFcn', num2str(Fs));

for i=1:9
    subsys_name =strcat('Ex2_E/Subsystem',num2str(i-1),'/');
    set_param(strcat(subsys_name,'L1'), 'l', num2str(M));
    set_param(strcat(subsys_name,'C1'), 'c', num2str(C));
    set_param(strcat(subsys_name,'R1'), 'R', num2str(R));
end

open_system("Ex2_E.slx", 'loadonly');
out = sim("Ex2_E.slx", dur);
disp("done")
vel0 = out.vel0.Data;
vel1 = out.vel1.Data;
vel2 = out.vel2.Data;
vel3 = out.vel3.Data;
vel4 = out.vel4.Data;
vel5 = out.vel5.Data;
vel6 = out.vel6.Data;

press = out.pressure.Data;
H0 = db(abs(fft(vel0)./fft(press)));
H1 = db(abs(fft(vel1)./fft(press)));
H2 = db(abs(fft(vel2)./fft(press)));
H3 = db(abs(fft(vel3)./fft(press)));
H4 = db(abs(fft(vel4)./fft(press)));
H5 = db(abs(fft(vel5)./fft(press)));
H6 = db(abs(fft(vel6)./fft(press)));

[val,f0] = findpeaks(H0);
f0 = 48000/144001 * f0;
%close all
figure('Renderer', 'painters', 'Position', [100 100 1000 500])
plot(f, H0, LineWidth=1.2);
xlim([0, 1600]); ylim([-120, -10])
hold on

plot(f, H1, LineWidth=1.2);
hold on

plot(f, H2, LineWidth=1.2);
hold on

plot(f, H3, LineWidth=1.2);
hold on

plot(f, H4, LineWidth=1.2);
hold on

plot(f, H5, LineWidth=1.2);
hold on

plot(f, H6, LineWidth=1.2);
hold on

grid minor
hold on
xline(f0(1), 'k--')
text(f0(1)*(1.01), -90, strcat('f_0=',num2str(f0(1))))

xline(f0(2), 'k--')
text(f0(2)*(1.01), -90, strcat('f_1=',num2str(f0(2))))

xline(f0(3), 'k--')
text(f0(3)*(1.01), -90, strcat('f_2=',num2str(f0(3))))

xline(f0(4), 'k--')
text(f0(4)*(1.01), -90, strcat('f_0=',num2str(f0(4))))

xline(f0(5), 'k--')
text(f0(5)*(1.01), -90, strcat('f_1=',num2str(f0(5))))

xline(f0(6), 'k--')
text(f0(6)*(1.01), -90, strcat('f_2=',num2str(f0(6))))

xline(f0(7), 'k--')
text(f0(7)*(1.01), -90, strcat('f_2=',num2str(f0(7))))

legend('H0', 'H1', 'H2','H3', 'H4', 'H5', 'H6')
xlabel('Freq [Hz]'); ylabel("|H| [dB]");
title("4x2 Tree, Frequency Response functions of the system")
%delete("./plots/Ex2_E_FRF.png");
%saveas(gcf, "./plots/Ex2_E_FRF.png");
