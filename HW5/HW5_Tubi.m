%% METODO DEI TUBI YEYE

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
deltaL = r1*0.85;
Lp = L+deltaL;
theta1 = atan(k0*x1);
theta2 = atan(k0*x2);

M = tan(k0*0.6*r2)*rho./k0./S2;

% Zin = 1i*rho*c./S2 .* sin(k0*Lp) .* sin(theta2) ./ sin(k0*Lp+theta2) + 1i*w0*M;
Zin = 1i*rho*c./S2 .* sin(k0*Lp) .* sin(theta1) ./ sin(k0*Lp+theta1) + 1i*w0*M;
Zindb = db(Zin);

figure()
plot(r1, Zindb, LineWidth=1.4)
grid minor
xlabel('r_f [m]'); ylabel('|Z_{in}| [dB]')

r1 = r1(Zindb == min(Zindb));
r2 = r1+L*tan(alpha);
deltaL = r1*0.6;
Lp = L+deltaL;

%% PUNTO B

f1 = 369.99;
w1 = 2*pi*f1;
k1 = w1/c;

x1 = r1/tan(alpha);
x2 = x1+L;
S1 = r1^2*pi;
S2 = r2^2*pi;
% delta = 0.6*r1;
theta1 = atan(k1*x1);
theta2 = atan(k1*x2);

M = tan(k1*0.6*r2)*rho./S2./k1;

D1 = linspace(0, L, 1000);
delta = D1 + deltaL.^2./(D1+2*deltaL);

Lpp = Lp-delta;

% Zin = 1i*w1*M + 1i*rho*c/S2 .* sin(k1*Lpp).*sin(theta2) ./ sin(k1*Lpp+theta2);
Zin = 1i*w1*M + 1i*rho*c/S2 .* sin(k1*Lpp).*sin(theta1) ./ sin(k1*Lpp+theta1);
Zindb = db(Zin);

figure()
plot(D1, Zindb, LineWidth=1.4)
grid minor
xlabel('D_1 [m]'); ylabel('|Z_{in}| [dB]')

hole3_coord = L-D1(Zindb == min(Zindb));
D1 = D1(Zindb == min(Zindb));
Lpp = Lp - (D1 + deltaL.^2./(D1+2*deltaL));

%% PUNTO C

f2 = 415.3;
w2 = 2*pi*f2;
k2 = w2/c;

D2 = linspace(0, Lpp, 1000);
delta2 = D2 - (D2*deltaL)./(D2+deltaL);

Lppp = Lpp-delta2;

% Zin = 1i*w2*M + 1i*rho*c/S2 .* sin(k2*Lppp).*sin(theta2) ./ sin(k2*Lppp+theta2);
Zin = 1i*w2*M + 1i*rho*c/S2 .* sin(k2*Lppp).*sin(theta1) ./ sin(k2*Lppp+theta1);
Zindb = db(Zin);

figure()
plot(D2, Zindb, LineWidth=1.4)
grid minor
xlabel('D_2 [m]'); ylabel('|Z_{in}| [dB]')

hole4_coord = Lpp-D2(Zindb == min(Zindb));
D2 = D2(Zindb == min(Zindb));
Lppp = Lpp - (D2 - (D2*deltaL)./(D2+deltaL));

%% shape plot

x = linspace(0, L, 1000);
shape = r2-x*tan(alpha);
hole1 = linspace(hole3_coord-r1, hole3_coord+r1, 10);
hole2 = linspace(hole4_coord-r1, hole4_coord+r1, 10);

r3 = (L-hole3_coord+x1)*tan(alpha);
r4 = (L-hole4_coord+x1)*tan(alpha);


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


