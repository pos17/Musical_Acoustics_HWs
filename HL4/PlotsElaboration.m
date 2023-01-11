%% HL4 Impedances plot and data manipulation
clc
clear all
close all
aspect_ratio = 0.7143; %fig aspect ratio

%% Zin Tube
clc
close all

T1 = readtable('Plots/Zin_Tube.csv','NumHeaderLines',5);  % skips the first three rows of data

figure('Renderer', 'painters', 'Position', [100 100 1200 600])
plot(T1{:,1}, T1{:,2},LineWidth=1.2)
[maxPks, maxLocs] = findpeaks(abs(T1{:,2})); % impedance maxima

[minPks, minLocs] = findpeaks(1-T1{:,2}); % impedance minima
hold on;
for ii=1:size(maxPks)
    xline(T1{maxLocs(ii),1},"--");
    text(T1{maxLocs(ii),1}-30, 50, strcat('f = ',num2str(T1{maxLocs(ii),1}), " Hz"),'FontSize',11);
    hold on
end

for ii=1:size(minPks)
    xline(T1{minLocs(ii),1},"--");
    text(T1{minLocs(ii),1}-30, 2, strcat('f = ',num2str(T1{minLocs(ii),1}), " Hz"),'FontSize',11);
    hold on
end

grid minor;
xlim([0,1270])
ylim([0,55])
xlabel("$ freq [Hz] $",'interpreter','latex','FontSize',12);
ylabel('$ \left | Z \right | [dB]$ ','interpreter','latex','FontSize',12)
title("Tube input impedance",'interpreter','latex','FontSize',15)

filename = "Zin_Tube";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);
%%
c=343;
n = 1:10;
Lt=1.37;
f_min = (n.*c)./(2*Lt);
f_max = ((2*n-1)*c)/(4*Lt);
%% Zin Tube Bell
T2 = readtable('Plots/Zin_TubeBell.csv','NumHeaderLines',5);  % skips the first three rows of data

figure('Renderer', 'painters', 'Position', [100 100 1200 600])
plot(T2{:,1}, T2{:,2},LineWidth=1.2)

%[maxPks, maxLocs] = findpeaks(abs(T{:,2})); % impedance maxima

%[minPks, minLocs] = findpeaks(1-T{:,2}); % impedance minima
%hold on;
%for ii=1:size(maxPks)
%    xline(T{maxLocs(ii),1},"--");
%    text(T{maxLocs(ii),1}-30, 6.5e-5, strcat('f = ',num2str(T{maxLocs(ii),1}), " Hz"),'FontSize',11);
%    hold on
%end

%for ii=1:size(minPks)
%    xline(T{minLocs(ii),1},"--");
%    text(T{minLocs(ii),1}-30, 0.5e-5, strcat('f = ',num2str(T{minLocs(ii),1}), " Hz"),'FontSize',11);
%    hold on
%end

grid minor;
xlim([0,1270])
xlabel("$ freq [Hz] $",'interpreter','latex','FontSize',12);
ylabel('$ \left | Z \right | [dB]$ ','interpreter','latex','FontSize',12)
title("Tube Bell input impedance",'interpreter','latex','FontSize',15)

filename = "Zin_TubeBell";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);

%% Zin Mouthpiece

figure('Renderer', 'painters', 'Position', [100 100 1200 600])
T3 = readtable('Plots/Zin_Mouthpiece.csv','NumHeaderLines',5);  % skips the first three rows of data
plot(T3{:,1}, T3{:,2},LineWidth=1.2)

[maxPks, maxLocs] = findpeaks(abs(T3{:,2})); % impedance maxima

%[minPks, minLocs] = findpeaks(1-T{:,2}); % impedance minima
hold on;
for ii=1:size(maxPks)
    xline(T3{maxLocs(ii),1},"--");
    text(T3{maxLocs(ii),1}+30, 62, strcat('f = ',num2str(T3{maxLocs(ii),1}), " Hz"),'FontSize',11);
    hold on
end

%for ii=1:size(minPks)
%    xline(T{minLocs(ii),1},"--");
%    text(T{minLocs(ii),1}-30, 0.5e-5, strcat('f = ',num2str(T{minLocs(ii),1}), " Hz"),'FontSize',11);
%    hold on
%end

grid minor;
xlim([0,1270])
xlabel("$ freq [Hz] $",'interpreter','latex','FontSize',12);
ylabel('$ \left | Z \right | [dB]$ ','interpreter','latex','FontSize',12)
title("Mouthpiece input impedance",'interpreter','latex','FontSize',15)

filename = "Zin_Mouthpiece";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);

%% Numerical Mouthpiece resonance frequency
rM = 0.009; %[m]
rT = 0.006; %[m]
rho = 1.225; %[kg/m^3]
c= 340;
l_c = 0.091;  %[m]
r_mean = ((rT/8)+rT)/2;
l_c_fixed = l_c + r_mean*(0.61+0.85);
S_c = (r_mean^2)*pi;

V = (4/3 * pi * rM^3 )
C = V/(rho*c^2)

L = (rho*l_c_fixed)/(S_c)

omega_0 = (L*C)^(-1/2)
f_0 = omega_0/(2*pi)

f_0_helmotz = (c/(2*pi))*sqrt(S_c/(V*l_c_fixed))


%% Zin Trumpet

figure('Renderer', 'painters', 'Position', [100 100 1200 600])
T4 = readtable('Plots/Zin_Trumpet.csv','NumHeaderLines',5);  % skips the first three rows of data
plot(T4{:,1}, T4{:,2},LineWidth=1.2)

%[maxPks, maxLocs] = findpeaks(abs(T4{:,2})); % impedance maxima

%[minPks, minLocs] = findpeaks(1-T{:,2}); % impedance minima
%hold on;
%for ii=1:size(maxPks)
%    xline(T3{maxLocs(ii),1},"--");
%    text(T3{maxLocs(ii),1}+30, 62, strcat('f = ',num2str(T3{maxLocs(ii),1}), " Hz"),'FontSize',11);
%    hold on
%end

%for ii=1:size(minPks)
%    xline(T{minLocs(ii),1},"--");
%    text(T{minLocs(ii),1}-30, 0.5e-5, strcat('f = ',num2str(T{minLocs(ii),1}), " Hz"),'FontSize',11);
%    hold on
%end

grid minor;
xlim([0,1270])
xlabel("$ freq [Hz] $",'interpreter','latex','FontSize',12);
ylabel('$ \left | Z \right | [dB]$ ','interpreter','latex','FontSize',12)
title("Trumpet input impedance",'interpreter','latex','FontSize',15)

filename = "Zin_Trumpet";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);

%% All graphs

figure('Renderer', 'painters', 'Position', [100 100 1200 600])
T4 = readtable('Plots/Zin_Trumpet.csv','NumHeaderLines',5);  % skips the first three rows of data
plot(T1{:,1}, T1{:,2},LineWidth=1.2)
hold on
plot(T2{:,1}, T2{:,2},LineWidth=1.2)
hold on
plot(T3{:,1}, T3{:,2},LineWidth=1.2)
hold on
plot(T4{:,1}, T4{:,2},LineWidth=1.2)

%[maxPks, maxLocs] = findpeaks(abs(T4{:,2})); % impedance maxima

%[minPks, minLocs] = findpeaks(1-T{:,2}); % impedance minima
%hold on;
%for ii=1:size(maxPks)
%    xline(T3{maxLocs(ii),1},"--");
%    text(T3{maxLocs(ii),1}+30, 62, strcat('f = ',num2str(T3{maxLocs(ii),1}), " Hz"),'FontSize',11);
%    hold on
%end

%for ii=1:size(minPks)
%    xline(T{minLocs(ii),1},"--");
%    text(T{minLocs(ii),1}-30, 0.5e-5, strcat('f = ',num2str(T{minLocs(ii),1}), " Hz"),'FontSize',11);
%    hold on
%end

grid minor;
legend(["Zin_{Tube}","Zin_{TubeBell}","Zin_{Mouthpiece}","Zin_{Trumpet}"])
xlim([0,1270])
xlabel("$ freq [Hz] $",'interpreter','latex','FontSize',12);
ylabel('$ \left | Z \right | [dB]$ ','interpreter','latex','FontSize',12)
title("Trumpet input impedance",'interpreter','latex','FontSize',15)

filename = "Zin_All";
delete([".\plots\"+filename+".png"]);
saveas(gcf, [".\plots\"+filename+".png"]);