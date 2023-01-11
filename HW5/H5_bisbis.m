clear; close all; clc;

L = 0.45; % m
alpha = deg2rad(0.75); % rad
c = 343; % m/s
f0 = 329.63; % Hz
w0 = 2*pi*f0;
rho = 1.225; % kg/m3
k0 = w0/c;

%% raius search

x1 = @(r) r/tan(alpha);
r_h = @(r) r+L*tan(alpha);
theta = @(r) atan(k0*x1(r));
deltaL = @(r) 0.6*r;
len = @(r) L+deltaL(r);
S_h = @(r) r_h(r)^2*pi;
S_f = @(r) r^2*pi;

M = @(r) deltaL(r)*rho/S_f(r);
Zm = @(r, w) 1i*w*M(r);
Zin = @(r) 1i*rho*c/S_h(r) * sin(k0*len(r)) * sin(theta(r)) / sin(k0*len(r)+theta(r));
Ztot = @(r) Zin(r)+Zm(r, w0);
Ztot_dB = @(r) db(Ztot(r));

r_F = fminsearch(Ztot_dB, 0.1);
r_H = r_h(r_F);

disp("R foot")
disp(r_F)
disp("R head")
disp(r_H)

xl = linspace(0, L, 100);
shape = r_H - xl*tan(alpha);

figure()
plot(xl, shape, 'k',LineWidth=1.2);
hold on
plot(xl, -shape, 'k', LineWidth=1.2)
hold on
yline(0, 'k--')
ylim([-0.1, 0.1])
grid minor
% xlabel("x [m]"); ylabel("|Z_{in}| [dB]")
title('Shape of the recorder resonator')

%%
rad = linspace(0, 0.1, 1000);
ZZZ = arrayfun(Ztot_dB, rad);

figure()
plot(rad, ZZZ, LineWidth=1.4)
hold on
xline(r_F, 'k--', LineWidth=1.4)
grid minor


%%
clc
f1 = 369.99;
w1 = 2*pi*f1;
k1 = w1/c;

x = linspace()

Zpar = @(Z1, Z2) Z1+Z2/(Z1*Z2);

% x20 = 

x20 = @(x) x+x1(r_F);
x10 = L+x1(r_F);
len_end = @(x) x1(r_F) + x + deltaL(r_F);
len_beg  = @(x) x10-x20(x);
r = @(x) (x+x1(r_F))*tan(alpha);
S = @(x)  r(x)^2*pi;
S_hole = S_f(r_F);
theta = @(x) atan(k1*x+k1*x1(r_F));



Zin = @(x) 1i*rho*c/S(x) * sin(k1*len_end(x)) * sin(theta(r)) / sin(k1*len_end(x)+theta(r));
Z_coneEnd = @(x) Zin(x) + Zm(r_F, w1);
Z_hole =@(x, k) 1i*rho*c/S_hole*tan(k*(r(x)+0.6*r_F));

Z_holeCone = @(x) Zpar(Z_coneEnd(x), Z_hole(x, k1));

num = @(x, k, length, ZL) rho*c/S(L) * ( 1i*ZL*sin(k*length-theta(x))/sin(theta(x)) + rho*c/S(x)*sin(k*length));
den = @(x, k, length, ZL) ZL*sin(k*length+theta(L)-theta(x))/sin(theta(x))/sin(theta(L)) - 1i*rho*c/S(x)*sin(k*length+theta(L))/sin(theta(L));

Ztot = @(x) num(x, k1, len_beg(x), Z_holeCone(x))/den(x, k1, len_beg(x), Z_holeCone(x));
Ztotdb = @(x) db(abs(Ztot(x)));

x_h1 = fminsearch(Ztot, x1(r_F));









