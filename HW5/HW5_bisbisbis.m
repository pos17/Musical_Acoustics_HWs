clear; close all; clc;

L = 0.45; % m
alpha = deg2rad(0.75); % rad
c = 343; % m/s
rho = 1.225; % kg/m3
f0 = 329.63; % Hz
w0 = 2*pi*f0;
k0 = w0/c;

%% PUNTO A

r1 = linspace(0, 0.06, 1000);
r2 = r1+L*tan(alpha);
x1 = r1/tan(alpha);
x2 = r2/tan(alpha);
S1 = r1.^2*pi;
S2 = r2.^2*pi;
Lp = L+r1*0.85;
theta1 = atan(k0*x1);
theta2 = atan(k0*x2);

M = tan(k0*0.6*r2)*rho./k0./S2;

Zin = 1i*rho*c./S2 .* sin(k0*Lp) .* sin(theta2) ./ sin(k0*Lp+theta2) + 1i*w0*M;
Zindb = db(Zin);

figure()
plot(r1, Zindb, LineWidth=1.4)
grid minor

r1 = r1(Zindb == min(Zindb));
r2 = r1+L*tan(alpha);


%% PUNTO B
% close all

f1 = 369.99;
% f1 = 392;
w1 = 2*pi*f1;
k1 = w1/c;

x1 = r1/tan(alpha);
x2 = r2/tan(alpha);

x3 = linspace(x1, x2, 1000);
x = linspace(0,L, 1000);

r3 = x3*tan(alpha);

S1 = r1^2*pi;
S2 = r2^2*pi;
S3 = r3.^2*pi;

% M = 0.04 * rho ./ S2;
M = tan(k1*0.6*r2)*rho./S2./k1;

theta1 = atan(k1*x1);
theta2 = atan(k1*x2);
theta3 = atan(k1*x3);

L23 = x2-x3;
L31 = x3-x1;
L31p = L31+0.85*r1;
% L31p = L31;
% l = r3+0.85*r1;
l = 0.6*r1;

Z_coneEnd = 1i*rho*c./S3 .* sin(k1*L31p) .* sin(theta3) ./ sin(k1*L31p + theta3);
Z_hole = 1i*rho*c/S1 .* tan(k1*l);
Zpar = @(Z1, Z2) (Z1.*Z2)./(Z1+Z2);

Z_ch = Zpar(Z_coneEnd, Z_hole);


num = 1i*Z_ch.*sin(k1*L23-theta3)./sin(theta3) + rho*c./S3 .* sin(k1*L23);
den = Z_ch.*sin(k1*L23+theta2-theta3)./sin(theta2)./sin(theta3) ...
    - 1i*rho*c./S3 .* sin(k1*L23+theta2)./sin(theta2);

Ztot = 1i*w1*M + rho*c/S2 * num./den;
Ztotdb = db(Ztot);

% x3 = x3(Ztotdb==min(Ztotdb));
x3 = x3(find(islocalmin(Ztotdb), 1));
hole3_coord = L-(x3-x1);

figure()
plot(x, flip(Ztotdb), LineWidth=1.4)
hold on
xline(hole3_coord, 'k--', LineWidth=1.4)

grid minor

%% PUNTO C

f2 = 415.3;
w2 = 2*pi*f2;
k2 = w2/c;

x = linspace(0, L, 1000);

x4 = linspace(x1,x2, 1000);

r4 = x4*tan(alpha);
r3 = x3*tan(alpha);

S3 = r3.^2*pi;
S4 = r4.^2*pi;

theta1 = atan(k2*x1);
theta2 = atan(k2*x2);
theta3 = atan(k2*x3);
theta4 = atan(k2*x4);

L43 = x4-x3;
L31 = x3-x1;
L24 = x2-x4;
L31p = L31+0.85*r1;
% l3 = r3+0.85*r1;
% l4 = r4+0.85*r1;
l3 = 0.85*r1;
l4 = l3;

M = tan(k2*0.6*r2)*rho./S2./k2;

Z_cone31 = 1i*rho*c./S3 .* sin(k2*L31p) .* sin(theta3) ./ sin(k2*L31p + theta3);
Z_hole3 = 1i*rho*c/S1 .* tan(k2*l3);
Zpar = @(Z1, Z2) (Z1.*Z2)./(Z1+Z2);

Z_3 = Zpar(Z_cone31, Z_hole3);


num = 1i*Z_3.*sin(k2*L43-theta3)./sin(theta3) + rho*c./S3 .* sin(k2*L43);
den = Z_3.*sin(k2*L43+theta4-theta3)./sin(theta3)./sin(theta4) ...
    - 1i*rho*c./S3 .* sin(k2*L43+theta4)./sin(theta4);

Z_cone43 = rho*c/S2 * num./den;
Z_hole4 = 1i*rho*c/S1 .* tan(k2*l4);

Z_4 = Zpar(Z_cone43, Z_hole4);


num = 1i*Z_4.*sin(k2*L24-theta4)./sin(theta4) + rho*c./S4 .* sin(k2*L24);
den = Z_4.*sin(k2*L24+theta2-theta4)./sin(theta2)./sin(theta4) ...
    - 1i*rho*c./S4 .* sin(k2*L24+theta2)./sin(theta2);

Zintot = j*w2*M + rho*c/S2 * num./ den;
Zintotdb = db(Zintot);

x4 = x4(find(islocalmin(Zintotdb), 1, 'first'));
r4 = x4*tan(alpha);
hole4_coord = L-(x4-x1);

figure()
plot(x, flip(Zintotdb), LineWidth=1.4)
hold on
xline(hole4_coord, 'k--', LineWidth=1.4)

grid minor

%% PLOT SHAPE

shape = r2-x*tan(alpha);
hole1 = linspace(hole3_coord-r1, hole3_coord+r1, 10);
hole2 = linspace(hole4_coord-r1, hole4_coord+r1, 10);


figure()
plot(x, shape, 'k', LineWidth=1.4);
hold on
plot(x, -shape, 'k', LineWidth=1.4);
hold on
yline(0, 'k--')
hold on 
xline(hole3_coord, 'b--');
hold on 
xline(hole4_coord, 'b--');
hold on
plot(hole1, ones(length(hole1), 1)*r3, LineWidth=5);
hold on
plot(hole2, ones(length(hole2), 1)*r4, LineWidth=5);
grid minor
ylim([-0.2 0.2])

 








