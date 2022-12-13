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
Npoints = flip(2:30);

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

%% PLOTTING ANALYTICAL VERSUS APPROXIMATED

close all 
figure('Renderer', 'painters', 'Position', [100 100 800 400])


[Zin,l] = eval_impedance(Npoints(1), L, a0, rho, c, k, Z);
plot(f, db(abs(Zin)),LineWidth=1.2,LineStyle="--",color="red")
hold on 
[Zin,l] = eval_impedance(Npoints(23), L, a0, rho, c, k, Z);
plot(f, db(abs(Zin)),LineWidth=1.2,LineStyle="--",color="green")
hold on
[Zin,l] = eval_impedance(Npoints(29), L, a0, rho, c, k, Z);
plot(f, db(abs(Zin)),LineWidth=1.2,color="black")
hold on 
 
plot(f, db(abs(Zin_anal)),LineStyle=":",LineWidth=1.2)
legend('Approximated 29 sections',  'Approximated 7 sections','Approximated 1 section','Analytical Impedance')
xlabel("Freq [Hz]"); ylabel("Impedance [dB]");
title("Analytical versus Approximated impedance")
grid minor
ylim([0,230])

filename = "Ex1AnVSAppr";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);


%% PLOTTING E1

%close all
figure('Renderer', 'painters', 'Position', [100 100 800 400])
loglog(deltas, e1, '-o',LineWidth=1.2)
hold on 
plot(deltas(1),e1(1),'o','MarkerSize',12,Color="red",LineWidth=1.4);
hold on 
plot(deltas(29),e1(29),'o','MarkerSize',12,Color="black",LineWidth=1.4);
hold on 
plot(deltas(23),e1(23),'o','MarkerSize',12,Color="green",LineWidth=1.4);
% xticks([deltas])
% xticklabels(deltas)
xticks(0.01:0.04:0.35)
xticklabels(0.01:0.04:0.35)
xlabel("\delta [m]"); ylabel("Mean Square Error");
title("Mean Square Error as function of \delta")

grid minor
filename = "Ex1A";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);

%% PLOTTING E2

figure('Renderer', 'painters', 'Position', [100 100 800 400])
plot(deltas, e2, '-o', LineWidth=1.2)
hold on 
plot(deltas(1),e2(1),'o','MarkerSize',12,Color="red",LineWidth=1.4);
hold on 
plot(deltas(29),e2(29),'o','MarkerSize',12,Color="black",LineWidth=1.4);
hold on 
plot(deltas(23),e2(23),'o','MarkerSize',12,Color="green",LineWidth=1.4);
xlabel("\delta [m]"); ylabel("freq error");
title("Frequency Error as function of \delta")
xticks(0.01:0.04:0.35)
xticklabels(0.01:0.04:0.35)
grid minor
hold on
filename = "Ex1B";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);

d0 = deltas(find(e2==min(e2), 1, 'last'));
N = Npoints(find(e2==min(e2), 1, 'last'));
N_cones = N-1;


%% Exponential Horn visualization
x = linspace(0,L,1000);

figure('Renderer', 'painters', 'Position', [100 100 800 400])
plot(x, a0*exp(m*x),'k','linewidth',1.5);
hold on;
plot(x, -a0*exp(m*x),'k','linewidth',1.5);
xlabel('Axial Distance [m]', 'Fontsize', 12)
ylabel('Radius[m]', 'Fontsize', 12)
ylim([-0.06 .06])
grid on;
title('Exponential Horn Profile', 'Fontsize', 12);
filename = "Ex1ExpHorn";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);
%% Exponential horn vs approximated 4 cones 

[coneX,coneY] = approximateShape(a0,m,Npoints(26),deltas(26));


figure('Renderer', 'painters', 'Position', [100 100 800 400])
plot(x, a0*exp(m*x),'r','linewidth',1.5);
hold on;

plot(coneX,coneY,'-o',LineWidth=1.2,Color="blue",LineStyle="--");
hold on
plot(x, -a0*exp(m*x),'r','linewidth',1.5);
hold on;
%plot(coneX, coneY,'o',Color'r','linewidth',1.5,LineStyle="--");
%hold on;
%plot(coneX, -coneY,'o',Color,'r','linewidth',1.5,LineStyle="--");
plot(coneX,-coneY,'-o',LineWidth=1.2,Color="blue",LineStyle="--");

xlabel('Axial Distance [m]', 'Fontsize', 12)
ylabel('Radius[m]', 'Fontsize', 12)
ylim([-0.06 .06])
grid on;

legend("Exponential Horn","4 sections approximated horn")
title('Exponential Horn Profile compared to approximated 4 sections ', 'Fontsize', 12);
filename = "Ex1ExpHorn4Appr";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);

%% %% Exponential horn vs three cases approximated 


[coneX_7,coneY_7] = approximateShape(a0,m,Npoints(23),deltas(23));
[coneX_29,coneY_29] = approximateShape(a0,m,Npoints(1),deltas(1));
[coneX_1,coneY_1] = approximateShape(a0,m,Npoints(29),deltas(29));

figure('Renderer', 'painters', 'Position', [100 100 800 400])
plot(x, a0*exp(m*x),'b','linewidth',1.5);
hold on;

plot(coneX_29,coneY_29,'-o',LineWidth=1.2,Color="red",LineStyle="--");
hold on
plot(coneX_7,coneY_7,'-o',LineWidth=1.2,Color="green",LineStyle="--");
hold on
plot(coneX_1,coneY_1,'-o',LineWidth=1.2,Color="black",LineStyle="--");
hold on


plot(x, -a0*exp(m*x),'b','linewidth',1.5);
plot(coneX_29,-coneY_29,'-o',LineWidth=1.2,Color="red",LineStyle="--");
hold on
plot(coneX_7,-coneY_7,'-o',LineWidth=1.2,Color="green",LineStyle="--");
hold on
plot(coneX_1,-coneY_1,'-o',LineWidth=1.2,Color="black",LineStyle="--");
xlabel('Axial Distance [m]', 'Fontsize', 12)
ylabel('Radius[m]', 'Fontsize', 12)
ylim([-0.06 .06])
grid on;

legend("Exponential Horn","29 sections approximated horn","7 sections approximated horn","1 section approximated horn")
title('Exponential Horn Profile compared to approximated 4 sections ', 'Fontsize', 12);
filename = "Ex1ExpHorn3Cases";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);
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
FRF_Plot('Zin', Zin2, f, max(f), h1,"Impedance considering radiation load")

filename = "Ex1CLoad";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);

% plot(f, db(abs(Zin2)), LineWidth=1.2)
% grid minor
h2 = figure('Renderer', 'painters', 'Position', [100 100 1000 500]);
FRF_Plot('Zin_{an}', Zin_anal, f, max(f), h2)
% plot(f, db(abs(Zin_anal)), LineWidth=1.2)
% grid minor
Z= zeros(1, length(f));
[Zin_ref,l] = eval_impedance(N, L, a0, rho, c, k, Z);

% computing error e_1 between Z_appr with and without load

arg = abs(Zin2-Zin_ref).^2;
int = trapz(w,arg);
e11_Load = 1/(w(end)-w(1))*int;

%% Comparison between impedance with and without load to analytical

figure('Renderer', 'painters', 'Position', [100 100 800 400])
subplot(2,1,1)
plot(f, db(abs(Zin2)),LineWidth=1.2)
hold on
plot(f, db(abs(Zin_ref)),LineWidth=1.2)
hold on
%plot(f, db(abs(Zin_anal)),LineWidth=1.2,LineStyle="--")
% hold on
grid minor
xlabel("freq [Hz]");
ylabel(['|','Z','| [dB]'])
title("Magnitude analytical vs approximated Impedance with load")
legend("Approximated Impedance","Analytical Impedance")
%legend("Impedance with load", "Impedance with no load")
subplot(2,1,2)
plot(f, angle(Zin2),LineWidth=1.2)
hold on
plot(f, angle(Zin_ref),LineWidth=1.2)
hold on
xlabel("freq [Hz]");
ylabel(['\angle','Z',' [rad]'])
%plot(f, angle(Zin_anal),LineWidth=1.2,LineStyle="--")
hold on
grid minor
title("Phase analytical vs approximated Impedance with load")
legend("Approximated Impedance","Analytical Impedance")
% legend("Impedance with load", "Impedance with no load")
filename = "Ex1CComp";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);
%% Comparison between analytical and approximated impedance with load

arg_AnApLoad = abs(Zin2-Zin_anal).^2;
int_AnApLoad = trapz(w,arg_AnApLoad);
e11_AnApLoad = 1/(w(end)-w(1))*int_AnApLoad;

figure('Renderer', 'painters', 'Position', [100 100 800 400])
subplot(2,1,1)
plot(f, db(abs(Zin2)),LineWidth=1.2)
hold on
plot(f, db(abs(Zin_anal)),LineWidth=1.2,LineStyle="--")
%hold on
%plot(f, abs(Zin_ref),LineWidth=1.2)
grid minor
xlabel("freq [Hz]");
ylabel(['|','Z','| [dB]'])
title("Magnitude analytical vs approximated Impedance with load")
legend("Approximated Impedance","Analytical Impedance")
%legend("Impedance with load", "Impedance with no load")
subplot(2,1,2)
plot(f, unwrap(angle(Zin2)),LineWidth=1.2)
hold on
% plot(f, angle(Zin_ref),LineWidth=1.2)
xlabel("freq [Hz]");
ylabel(['\angle','Z',' [rad]'])
plot(f, angle(Zin_anal),LineWidth=1.2,LineStyle="--")
%hold on
grid minor
title("Phase analytical vs approximated Impedance with load")
legend("Approximated Impedance","Analytical Impedance")
% legend("Impedance with load", "Impedance with no load")
filename = "Ex1CCompLoad";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);
%% Comparison between analytical and approximated impedance with no load

ZL = 0;
num_ref_anal = ZL.*cos(b*L+theta) + 1i*rho*c/S2_anal*sin(b*L);
den_ref_anal = 1i*ZL.*sin(b*L) + rho*c/S2_anal * cos(b*L-theta);
Zin_ref_anal = rho*c/S1_anal*(num_ref_anal./den_ref_anal);

arg_AnApnoLoad = abs(Zin_ref-Zin_ref_anal).^2;
int_AnApnoLoad = trapz(w,arg_AnApnoLoad);
e11_AnApnoLoad = 1/(w(end)-w(1))*int_AnApnoLoad;

figure('Renderer', 'painters', 'Position', [100 100 800 400])

subplot(2,1,1)
plot(f, db(abs(Zin_ref)),LineWidth=1.2)
hold on
plot(f, db(abs(Zin_ref_anal)),LineWidth=1.2,LineStyle="--")
grid minor
xlabel("freq [Hz]");
ylabel(['|','Z','| [dB]'])
title("Magnitude analytical vs approximated Impedance without load")
legend("Approximated Impedance","Analytical Impedance")
%legend("Impedance with load", "Impedance with no load")
subplot(2,1,2)
plot(f, angle(Zin_ref),LineWidth=1.2)
hold on
plot(f, angle(Zin_ref_anal),LineWidth=1.2,LineStyle="--")
xlabel("freq [Hz]");
ylabel(['\angle','Z',' [rad]']);
grid minor
title("Phase analytical vs approximated Impedance without load")
legend("Approximated Impedance","Analytical Impedance")
filename = "Ex1CCompNoLoad";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);

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
f_res10 = f(locs(1:10));

%% Computing compound horn plot
pipeX = [0,0.85];
pipeY = [0.008,0.008];
coneX_7_Compound = coneX_7 +0.85;
figure('Renderer', 'painters', 'Position', [100 100 800 400])

plot(pipeX,pipeY,'k',LineWidth=1.2)
hold on
plot(pipeX,-pipeY,'k',LineWidth=1.2)
hold on 
plot(coneX_7_Compound,coneY_7,'-o',LineWidth=1.2,Color="black")
hold on 
plot(coneX_7_Compound,-coneY_7,'-o',LineWidth=1.2,Color="black")
xlabel('Axial Distance [m]', 'Fontsize', 12)
ylabel('Radius[m]', 'Fontsize', 12)
ylim([-0.06 .06])
grid on;
title('Compound Horn Profile', 'Fontsize', 12);
filename = "Ex2DCompoundHorn";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);


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



