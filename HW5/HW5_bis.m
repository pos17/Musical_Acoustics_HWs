clear; close all; clc;

L = 0.45; % m
alpha = deg2rad(0.75); % rad
c = 343; % m/s
f0 = 329.63; % Hz
w0 = 2*pi*f0;
rho = 1.225; % kg/m3
k0 = w0/c;

%% DIAMETERS
clc
% r_h = linspace(L*tan(alpha), 0.1, 1000);
r_f = linspace(0, 0.06, 10000); % foot
r_h = r_f+L*tan(alpha); % head 
x1 = r_f/tan(alpha); 
theta1 = atan(k0*x1);
len = L + 0.6*r_f;
S_h = r_h.^2 * pi;
S_f = r_f.^2 * pi;

Zin = 1i*rho*c./S_h .* sin(k0*len).*sin(theta1) ./ sin(k0*len + theta1);
deltaL = 0.6*r_f; % come da fletching
% deltaL = 40e-3; % se vuoi
M = deltaL * rho ./ S_f; % come da fletching, contrario di fede
Zm = 1i*w0*M;
Ztot = Zin+Zm;
Ztotdb = db(abs(Ztot));

close all
figure()
% plot(r_h, db(abs(Zin)), LineWidth=1.2)
% hold on
plot(r_f, db(abs(Ztot)), LineWidth=1.2)
grid minor
xlabel("r_{foot} [m]"); ylabel("|Z_{in}| [dB]")
title('Impedance as function of r_f')

r_h_found = r_h(Ztotdb==min(Ztotdb));
% r_f_found = r_f(islocalmin(db(abs(Ztot))));
% r_f = r_f_found(3)
r_f = r_f(Ztotdb==min(Ztotdb))
r_h = r_f+L*tan(alpha)

% x1 = r_f_found/tan(alpha)


%% SHAPE PLOTTING

x = linspace(0, L, 100);
shape = r_h - x*tan(alpha);

figure()
plot(x, shape, 'k',LineWidth=1.2);
hold on
plot(x, -shape, 'k', LineWidth=1.2)
hold on
yline(0, 'k--')
ylim([-0.1, 0.1])
grid minor
% xlabel("x [m]"); ylabel("|Z_{in}| [dB]")
title('Shape of the recorder resonator')

%% IMPEDANCE PLOT 

f = linspace(200, 800, 1000);
w = 2*pi*f;
k = w./c;


x1 = r_f/tan(alpha); 
theta1 = atan(k0*x1);
len = L + 0.6*r_f;
S_h = r_h.^2 * pi;
S_f = r_f.^2 * pi;

Zin = 1i*rho*c./S_h .* sin(k*len).*sin(theta1) ./ sin(k*len + theta1);
deltaL = 40e-3;
deltaL = 0.6*r_f; % come da fletching
M = deltaL * rho ./ S_f; % come da fletching, contrario di fede
Zm = 1i*w*M;
Ztot = Zin+Zm;


figure()
plot(f,db(abs(Ztot)), LineWidth=1.2)
grid minor
xlabel("freq [Hz]"); ylabel("|Z_{in}| [dB]")
title('Impedance as function of freq')

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







