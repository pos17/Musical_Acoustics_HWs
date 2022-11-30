close all; clear; clc;

% INITIAL PARAMETERS
a0 = 0.008; % m
m = 4.2; % m^-1
L = 0.35; % m

c = 343; % m/s
rho = 1.21; % kg/m3


%freq axis
f = linspace(1, 2000, 2000);
w = 2*pi*f;
k = w/c;

%% ANALYTICAL IMPEDANCE
b = sqrt(k.^2-m^2);
theta = atan(m./b);
S2_anal = (a0*exp(m*L))^2*pi;
S1_anal = a0^2*pi;

num_anal = 1i*rho*c/S2_anal*sin(b*L);
den_anal = rho*c/S2_anal * cos(b*L-theta);
Zin_anal = rho*c/S1_anal*(num_anal./den_anal);

%% ERRORS
% Npoints = [0.1, 0.2, 0.5, 0.8, 1, 2, 5, 8, ... 
%     10, 20, 50, 80, 100, 200, 500, 800, ...
%     1000, 2000, 5000]*100;

% Npoints = flip(round(logspace(0, 2, 20)));
Npoints = flip(1:30);

e1 = zeros(1, length(Npoints));
e2 = zeros(1, length(Npoints));
deltas = zeros(1, length(Npoints));
Z= zeros(1, length(f));
for jj = 1:length(e1)
    [Zin,l] = eval_impedance(Npoints(jj), L, a0, rho, c, k, Z);
    arg = abs(Zin-Zin_anal).^2;
    int = trapz(w,arg);
    e1(jj) = 1/(w(end)-w(1))*int;
    
    [pks1, locs1] = findpeaks(abs(Zin));
    [pks2, locs2] = findpeaks(abs(Zin_anal));

    for ii=1:length(locs1)
         e2(jj)=e2(jj)+min(abs(w(locs1(ii))-w(locs2(ii))));
    end
    
    deltas(jj) = l;
end

%% PLOTTING E1

close all
figure('Renderer', 'painters', 'Position', [100 100 800 400])
semilogy(deltas, e1, LineWidth=1.2)
xlabel("\delta [m]"); ylabel("Mean Square Error");
title("Mean Square Error as function of \delta")
grid minor
% hold on
% xline(d0, 'k--');
% text(d0*1.5,min(e1), ["\delta="+num2str(d0)] )
filename = "Ex1A";
% delete([".\plots\"+filename+".png"]);
% saveas(gcf, [".\plots\"+filename+".png"]);

%% PLOTTING E2

figure('Renderer', 'painters', 'Position', [100 100 800 400])
plot(deltas, e2, LineWidth=1.2)
xlabel("\delta [m]"); ylabel("freq error");
title("Frequency Error as function of \delta")
grid minor
hold on
filename = "Ex1B";
% delete([".\plots\"+filename+".png"]);
% saveas(gcf, [".\plots\"+filename+".png"]);

d0 = deltas(find(e2==min(e2), 1, 'last'));
N = Npoints(find(e2==min(e2), 1, 'last'));
N_cones = N-1;

%% POINT C
e22 = 0;
% Zin = eval_impedance(N, L, a0, rho, c, k, zeros(1,length(k)));
[Zin2, ZL] = eval_impedance2(N, L, a0, rho, c, k);

num_anal = ZL.*cos(b*L+theta) + 1i*rho*c/S2_anal*sin(b*L);
den_anal = 1i*ZL.*sin(b*L) + rho*c/S2_anal * cos(b*L-theta);
Zin_anal = rho*c/S1_anal*(num_anal./den_anal);

arg = abs(Zin_anal-Zin2).^2;
int = trapz(w,arg);
e11 = 1/(w(end)-w(1))*int;

[pks1, locs1] = findpeaks(abs(Zin2));
[pks2, locs2] = findpeaks(abs(Zin_anal));

figure
plot(f, db(abs(Zin2)), LineWidth=1.2)
grid minor
figure
plot(f, db(abs(Zin_anal)), LineWidth=1.2)
grid minor

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% P A R T  2
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% POINT D

L1 = 0.85; % m
Z0 = rho*c/S1_anal;
num = Zin2.*cos(k*L1) + 1i*Z0.*sin(k*L1);
den = 1i*Zin2.*sin(k*L1) + Z0.*cos(k*L1);
Z_comp = Z0.*num./den;

close all;
figure;
plot(f, db(abs(Z_comp)));
grid minor

[pks, locs] = findpeaks(abs(Z_comp));
f_res10 = f(locs(1:10));

%% POINT E


