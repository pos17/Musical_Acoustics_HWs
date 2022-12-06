close all; clear; clc;

% INITIAL PARAMETERS
a0 = 0.008; % m
m = 4.2; % m^-1
L = 0.35; % m

c = 343; % m/s
rho = 1.21; % kg/m3


%freq axis
% f = linspace(1, 2000, 2000);
f = 1:2000;
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
    
    [pks10, locs10] = findpeaks(abs(Zin));
    [pks20, locs20] = findpeaks(abs(Zin_anal));

    for ii=1:length(locs10)
         e2(jj)=e2(jj)+min(abs(w(locs10(ii))-w(locs20(ii))));
    end
    
    deltas(jj) = l;
end

%% PLOTTING E1

close all
figure('Renderer', 'painters', 'Position', [100 100 800 400])
loglog(deltas, e1, '-o',LineWidth=1.2)
xlabel("\delta [m]"); ylabel("Mean Square Error");
title("Mean Square Error as function of \delta")
grid minor
filename = "Ex1A";
% delete([".\plots\"+filename+".png"]);
% saveas(gcf, [".\plots\"+filename+".png"]);

%% PLOTTING E2

figure('Renderer', 'painters', 'Position', [100 100 800 400])
plot(deltas, e2, '-o', LineWidth=1.2)
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

[pks11, locs11] = findpeaks(abs(Zin2));
[pks21, locs21] = findpeaks(abs(Zin_anal));

h1 = figure('Renderer', 'painters', 'Position', [100 100 1000 500]);
FRF_Plot('Zin', Zin2, f, max(f), h1)
% plot(f, db(abs(Zin2)), LineWidth=1.2)
% grid minor
h2 = figure('Renderer', 'painters', 'Position', [100 100 1000 500]);
FRF_Plot('Zin_{an}', Zin_anal, f, max(f), h2)
% plot(f, db(abs(Zin_anal)), LineWidth=1.2)
% grid minor

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

% close all
jj = figure;
% plot(f, db(abs(Z_comp)), line);
FRF_Plot('Z', Z_comp, f, 2000, jj)


[pks, locs] = findpeaks(abs(Z_comp));
f_res10 = f(locs(1:10))';

%% POINT E - Inharmonicity

% close all
inh_exp = f(locs10);
inh_exp2 = f(locs11);
inh_exp = inh_exp./inh_exp(1)./(1:length(inh_exp));
inh_exp2 = inh_exp2./inh_exp2(1)./(1:length(inh_exp2));

% f_res10 = f_res10
inh_comp = f_res10./f_res10(1)./(1:length(f_res10));

disp("Horn inharmonicity")
disp(inh_exp)
disp("Horn w/ ext impedance inharmonicity")
disp(inh_exp)
disp("Compound Horn inharmonicity")
disp(inh_comp);

figure('Renderer', 'painters', 'Position', [100 100 1000 500])
plot(inh_exp, '-o', LineWidth=1.2)
hold on
plot(inh_exp2, '-o', LineWidth=1.2)
hold on 
plot(inh_comp, '-o', LineWidth=1.2)
hold on 
plot(1:10, ones(1,10), '-o' ,LineWidth=0.8)
grid
legend(["horn", "horn w/ load", "compound horn", "perfect harmonicity"], 'Position', [0.717914508303892,0.536260275699164,0.172624437129336,0.142886516294981])
xlabel('# of harmonic'); ylabel('harmonicity')



