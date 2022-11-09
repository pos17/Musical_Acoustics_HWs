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

%% EX2 - electrical

M =  rho*l/S;
C = V0/(rho*c^2);
R = rho*c/S;

f = linspace(0, Fs, N);
w = 2*pi*f;

Z = R+j*w*M+(j*w*C).^(-1);

set_param('Ex2', 'PreLoadFcn', num2str(Fs));


%%

for ii=1:7
    set_param(['Ex2/L' num2str(ii)], 'l', num2str(M));
    set_param(['Ex2/C' num2str(ii)], 'c', num2str(C));
    set_param(['Ex2/R' num2str(ii)], 'R', num2str(R));
end

open_system("HL2\Ex2.slx", 'loadonly');
out = sim("HL2\Ex2.slx", dur);

%% plotting
vel0 = out.vel0.Data;
vel1 = out.vel1.Data;
vel2 = out.vel2.Data;
press = out.pressure.Data;

H0 = db(abs(fft(vel0)./fft(press)));
H1 = db(abs(fft(vel1)./fft(press)));
H2 = db(abs(fft(vel2)./fft(press)));

f0 = f(find(H0==max(H0),1));
f1 = f(find(H1==max(H1),1));
f2 = f(find(H2==max(H2),1));

close all
figure;
plot(f, H0, LineWidth=1.2);
xlim([0, f0*2]); ylim([-100, 0])
hold on
plot(f, H1, LineWidth=1.2);
% xlim([0, Fs/2]);
hold on
plot(f, H2, LineWidth=1.2);
% xlim([0, Fs/2]);
grid minor
hold on
xline(f0, 'k--')
legend('H0', 'H1', 'H2')

%%
close all

figure;
for nn = 1:5
    Ztot = 0;
    N=nn;
    K=2;
    for ii=0:N
        Ztot = Ztot + Z/(K)^ii;
    end
    
    mob = 1./Ztot;
    plot(f, db(abs(mob)), LineWidth=1.2);
    xlim([0, f0*2]); ylim([-100, 0])
    hold on
end
legend('N=1','N=2','N=3','N=4','N=5')
grid minor


figure;
for kk = 2:6
    Ztot = 0;
    N=3;
    K=kk;
    for ii=0:N
        Ztot = Ztot + Z/(K)^ii;
    end
    
    mob = 1./Ztot;
    plot(f, db(abs(mob)), LineWidth=1.2);
    xlim([0, f0*2]); ylim([-100, 0])
    hold on
end
legend('K=1','K=2','K=3','K=4','K=5')
grid minor
