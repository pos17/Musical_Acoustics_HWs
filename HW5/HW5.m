clear; close all; clc;

L = 0.45; % m
alpha = deg2rad(0.75); % grad
c = 343; % m/s
f0 = 329.63; % Hz
w0 = 2*pi*f0;
rho = 1.21; % kg/m3
k0 = w0/c;

%% DIAM 

%% 1) DIAMETERS
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x1 = r1/tan(alpha);                                                     %  
% S1 = pi*r1^2;                                                           %  
% theta1 = atan(k*x1);                                                    %  
%                                                                         %  
% Zin = abs(1i*(rho*c/S1) * sin(k*L)*sin(theta1)/sin(k*L+theta1));        %
%                                                                         %  
% we want to the denominator 'sin(k*L+theta1))' to be zero                %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x20 = @(r) r/tan(alpha);
x10 = @(r) x20(r)+L;
R = @(r) r+L*tan(alpha); % aperture diameter

theta = @(x) atan(k0*x);
S = @(r) r^2*pi;


% Sp = @(r) pi*R(r)^2;
% Ss = @(r) Sp(r)/(1+cos(alpha));
% Z_L0 = @(r) 0.25*w0^2*rho/pi/c + 1i*0.61*rho*w0/pi/R(r);
% Z_rad = @(r) Z_L0(r)*Sp(r)/Ss(r);

Sp = @(r) pi*r^2;
Ss = @(r) Sp(r)/(1+cos(alpha));
Z_L0 = @(r) 0.25*w0^2*rho/pi/c + 1i*0.61*rho*w0/pi/r;
Z_rad = @(r) Z_L0(r)*Sp(r)/Ss(r);

num0 = @(r, k, ZL) rho*c/S(R(r)) * 1i*ZL*( sin(k*L-theta(x20(r)))/sin(theta(x20(r))) ) + rho*c/S(r) * sin(k*L);
den0 = @(r, k, ZL) ZL*( sin(k*L+theta(x10(r))-theta(x20(r)))/( sin(theta(x10(r)))*sin(theta(x20(r))) ) ) -1i*rho*c/S(r)*sin(k*L+theta(x10(r)))/sin(theta(x10(r)));
Y_cone = @(r) db(abs(den0(r, k0, Z_rad(r))/num0(r, k0, Z_rad(r))));

% mouth diameter
r2 = fminsearch(Y_cone, 0);
% throat diameter
r1 = R(r2);

close all
figure();
rad = linspace(0, 0.01, 1000);
ciao = arrayfun(Y_cone, rad);
plot(rad, (ciao), LineWidth=1.5)
hold on
xline(r2, 'r--', LineWidth=1)
xline(rad(ciao==min(ciao)), 'k--')

diam1 = 2*r1;
diam2 = 2*r2;

% disp("Throat diameter: ")
% disp(["   "+num2str(diam1*1e3)+" mm"])
% disp("Mouth diameter: ")
% disp(["   "+num2str(diam2*1e3)+" mm"])
disp("Throat rad: ")
disp(["   "+num2str(r1)+" m"])
disp("Mouth rad: ")
disp(["   "+num2str(r2)+" m"])
%% 2) FIRST HOLE

f1 = 369.29; % Hz
w1 = 2*pi*f1;
k1 = w1/c;
x1 = x10(r1);
x_end = x20(r1);

Z_L_rad = Z_rad(r1);
Z_par = @(Z1, Z2) Z1*Z2/(Z1+Z2);
theta = @(x) atan(k1*x);
L_fin = @(x) x1+L-x;
rd = @(x) (x)*tan(alpha);
l = @(x) rd(x)+0.6*diam2/2;
S = @(x) rd(x)^2*pi;

Z_h = @(x, k) 1i*rho*c/S(L+x1)*tan(k*l(x));
% Z_coneEnd = @(x, k) 1i*rho*c/S(x)*sin(k*L_fin(x))*sin(theta(x))/sin(k*L_fin(x)+theta(x));
num01 = @(r, k, ZL, len) rho*c/S(r) * 1i*ZL*( sin(k*len-theta(x20(r)))/sin(theta(x20(r))) ) + rho*c/S(R(r1)) * sin(k*len);
den01 = @(r, k, ZL, len) ZL*( sin(k*len+theta(x10(r))-theta(x20(r)))/( sin(theta(x10(r)))*sin(theta(x20(r))) ) ) - 1i*rho*c/S(R(r1))*sin(k*len+theta(x10(r)))/sin(theta(x10(r)));
Z_coneEnd = @(x, k) num01(rd(x), k0, Z_L_rad, x_end-x)/den01(rd(x), k0, Z_L_rad, x_end-x);
Z_coneHole = @(x, k) Z_par(Z_h(x, k), Z_coneEnd(x, k));

num1 = @(x, k, ZL) rho*c/S(x1) * 1i*ZL*( sin(k*(x-x1)-theta(x))/sin(theta(x)) ) + rho*c/S(x) * sin(k*(x-x1));
den1 = @(x, k, ZL) ZL*( sin(k*(x-x1)+theta(x1)-theta(x))/( sin(theta(x1))*sin(theta(x)) ) ) -1i*rho*c/S(x)*sin(k*(x-x1)+theta(x1))/sin(theta(x1));

Y_tot1 = @(x) den1(x, k1, Z_coneHole(x, k1))/num1(x, k1, Z_coneHole(x, k1));
Y_tot1_ABS = @(x) abs(Y_tot1(x));
Z_tot1 = @(x) 1/Y_tot1(x);

x2 = fminsearch(Y_tot1_ABS, 0);
x_H1 = x2-x1;

disp("First hole position: ")
disp(["   "+num2str(x_H1)+" m"])

%% 2) SECOND HOLE
f2 = 415.3;
w2 = 2*pi*f2;
k2 = w2/c;

num2 = @(x, k, ZL) rho*c/S(x) * 1i*ZL*( sin(k*(x2-x)-theta(x2))/sin(theta(x2)) ) + rho*c/S(x) * sin(k*(x2-x));
den2 = @(x, k, ZL) ZL*( sin(k*(x2-x)+theta(x)-theta(x2))/( sin(theta(x))*sin(theta(x2)) ) ) - 1i*rho*c/S(x2)*sin(k*(x2-x)+theta(x))/sin(theta(x));

Zin_cone2 = @(x) num2(x, k2, Z_coneHole(x2, k2))/den2(x, k2, Z_coneHole(x2, k2));
Z_coneHole2 = @(x) Z_par(Z_h(x, k2), Zin_cone2(x));

Y_tot2 = @(x) den1(x, k2, Z_coneHole2(x))/num1(x, k2, Z_coneHole2(x));
Y_tot2_ABS = @(x) abs(Y_tot2(x));
Z_tot2 = @(x) 1/Y_tot2(x);

x3 = fminsearch(Y_tot2_ABS, 0);
x_H2 = x3-x1;


disp("Second hole position: ")
disp(["   "+num2str(x_H2)+" m"])

%% PLOTTING HOLES + ADMITTANCE
close all
figure()
x = linspace(0, L, 10000);
Y_T1 = arrayfun(Y_tot1, x+x1);
plot(x, pow2db(abs(Y_T1)), LineWidth=1.5);
hold on
grid minor
xline(x_H1, 'k--')
xlabel('x [m]'); ylabel('|Y(x)| [dB]')
title('Admittance for first hole', 'f_0 = 369.29 Hz')

figure()
Y_T2 = arrayfun(Y_tot2, x+x1);
plot(x, pow2db(abs(Y_T2)), LineWidth=1.5);
hold on
grid minor
xline(x_H2, 'k--')
xlabel('x [m]'); ylabel('|Y(x)| [dB]')
title('Admittance for second hole', 'f_0 = 415.3 Hz')

%% SHAPE PLOTTING

x = linspace(0, L, 1000);
shape =  diam2/2 - tan(alpha)*x;

close all
figure()
plot(x, shape, 'k', LineWidth=1.5);
hold on
plot(x, -shape, 'k', LineWidth=1.5);
hold on
yline(0, 'k--', LineWidth=1.1)
hold on 
xline(x_H1, 'k--')
hold on 
xline(x_H2, 'k--')
grid minor
ylim([-0.02, 0.02])

%% 4) FLUE CHANNEL THICKNESS

fc = 1.7e3; % Hz, centroid
delta_P = 55; % Pa

Uj = sqrt(2*delta_P/rho);
h = 0.3*Uj/fc;

nu = 1.5e-5; % m2/s
Re = Uj*h/nu;

disp("Flue thickness:")
disp(["   "+num2str(h*1e3)+" mm"])
disp("Reynolds number:")
disp(["   "+num2str(Re)])
%% 5) BUONDARY LAYER

l_flue = 25e-3; % m 
d = sqrt(nu*l_flue/Uj);

disp("Boundary layer thickness: ")
disp(["   "+num2str(d*1e3)+" mm"])







