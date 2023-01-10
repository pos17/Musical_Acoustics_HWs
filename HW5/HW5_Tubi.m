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
deltaL = r1*0.6;
Lp = L+deltaL;
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

D = linspace(0, L, 1000);
delta = D + deltaL.^2./(D+2*deltaL);

Lpp = Lp-delta;

Zin = 1i*w1*M + 1i*rho*c/S2 .* sin(k1*Lpp).*sin(theta2) ./ sin(k1*Lpp+theta2);
Zindb = db(Zin);

figure()
plot(D, Zindb, LineWidth=1.4)
grid minor

hole3_coord = L-D(Zindb == min(Zindb));
D = D(Zindb == min(Zindb));
Lpp = Lp - (D + deltaL.^2./(D+2*deltaL));

%% PUNTO C

f2 = 415.3;
w2 = 2*pi*f2;
k2 = w2/c;

D2 = linspace(0, Lpp, 1000);
delta2 = D2 - (D2*deltaL)./(D2+deltaL);

Lppp = Lpp-delta2;

Zin = 1i*w2*M + 1i*rho*c/S2 .* sin(k2*Lppp).*sin(theta2) ./ sin(k2*Lppp+theta2);
Zindb = db(Zin);

figure()
plot(D2, Zindb, LineWidth=1.4)
grid minor

hole4_coord = Lpp-D2(Zindb == min(Zindb));
D2 = D2(Zindb == min(Zindb));
Lppp = Lpp - (D + deltaL.^2./(D+2*deltaL));


