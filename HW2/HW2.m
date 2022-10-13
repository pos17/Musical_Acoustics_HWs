%%%%%%% HOMEWORK 2 %%%%%%%%

clc; close all; clear;

addpath(genpath('HW2'))
T = 10; % N/m
a = 0.15; % m
sigma = 0.07; % kg/m^2

%% PART 1 - a

c = sqrt(T/sigma);

%% PART 1-b

J.values = [
    2.4048, 3.8317, 5.1356, ...
    5.5201, 6.3802, 7.0156, ...
    7.5883, 8.4172, 8.6537, ...
    8.7715, 9.7610, 9.9361, ...
    10.1735, 11.0647, 11.0864, ...
    11.6198, 11.7915, 12.2251
    ];

J.idx = {
    [0,1], [1,1], [2,1], ...
    [0,2], [3,1], [1,2], ...
    [4,1], [2,2], [0,3], ...
    [5,1], [3,2], [6,1], ...
    [1,3], [4,2], [7,1], ...
    [2,3], [0,4], [8,1]
    };

freqs = zeros(18, 1);
for i=1:length(J.values)
    freqs(i) = J.values(i)/(2*pi*a)*c;
    disp(['freq: ', num2str(freqs(i)), ', mn: ', num2str(cell2mat(J.idx(i)))] )
end


theta = linspace(0, 2*pi, 101);
r = linspace(0, a, 101);

close all
figure('Renderer', 'painters', 'Position', [10 10 1200 800])
for i=1:6
    subplot(2,3, i)
    idx = cell2mat(J.idx(i));
    m = idx(1);
    kn = freqs(i)*2*pi/c;
    Jm = besselj(m, kn*r);
    z = exp(1j*m*theta)'*Jm;
    polarplot3d(real(z'), 'RadialRange', [0, a], 'TickSpacing', 0);
    view(2);
    colorbar
    grid minor
    title(['f_{ ', num2str(idx), '}= ', num2str(freqs(i)), ' Hz'], Interpreter="tex")
end
sgtitle("First six modes of resonance")


%% PART 1- c

f = linspace(1,200, 4000);
omega = 2*pi*f;
Q = 25;
modes = 18;
graphmodes = 5;

close all
% figure
H = zeros(18, length(f));

figure('Renderer', 'painters', 'Position', [10 10 1200 800])

for i = 1:modes
    for j = 1:modes
        omega0i = (2*pi*freqs(i));
        H(i, :) = 1./(-omega.^2 + 1j*omega*omega0i/Q + omega0i^2);
        if(j==i && (i==1 || i==2 || i==3 || i==modes))
            idx = i;
            if(i==modes) 
                idx=graphmodes;
            end
            subplot(graphmodes, graphmodes, sub2ind([graphmodes, graphmodes], idx, idx));
            plot(f, abs(H(i,:)));
            grid minor
            title(['mode', num2str(i)])
        end
    end
end
subplot(graphmodes, graphmodes, sub2ind([graphmodes, graphmodes], graphmodes, 1));
plot(f, zeros(1, length(f)));
grid minor
subplot(graphmodes, graphmodes, sub2ind([graphmodes, graphmodes], 1, graphmodes));
plot(f, zeros(1, length(f)));
grid minor


close all

x = zeros(2, modes);
% phi = [15, 195];
radius = 0.075;

for nn = 1:2
    for i=1:modes
        idx = cell2mat(J.idx(i));
        m = idx(1);
        kn = freqs(i)*2*pi/c;
        Jm = besselj(m, kn*r);
        Jm = Jm(r==radius);
        if m~=0
            angle = pi/m + (nn-1)*pi;
            x(nn, i) = exp(1j*m*(angle*2*pi/360))*Jm;
        else
            angle = (nn-1)*pi;
            x(nn, i) = exp(1j*m*(angle*2*pi/360))*Jm;
        end
    end
end

H21 = zeros(1,length(f));
mob = zeros(1,length(f));

for ii = 1:length(f)
    zer = zeros(1, 200);
    H_tmp = diag(H(:,ii));
    H21(ii) = x(1,:)*H_tmp*x(2,:)';
    mob(ii) = 1j*2*pi*f(ii)*H21(ii);
end

figure
plot(f, abs(H21), LineWidth=1.2, LineWidth=1.2);
grid minor
title('Module of FRF_{2,1}', Interpreter='tex');



len = 2*length(f);
t = linspace(0, 50, len);
force = 0.1*exp(-(t-0.03).^2/0.01^2);
force_omega = fft(force);
force_omega = force_omega(1:len/2);

% close all
figure
plot(f, abs(force_omega), LineWidth=1.2);
grid minor
title('F(\omega)', Interpreter='tex')

% displacement = force_omega.*mob;
displacement = force_omega.*H21;

% close all
figure
plot(f, abs(displacement), LineWidth=1.2);
grid minor
title('Z(\omega)', Interpreter='tex')


dsp = ifft(displacement, len);

% close all
figure
plot(t, real(dsp), LineWidth=1.2);
% hold on
% plot(t, -force*1e-5, LineWidth=1.2);
grid minor
xlim([0, 4])
title('z(t)')




%% Part 2) 
E = 69e+9;  %[Pa]
rho = 2700; %[kg/m^3]
v= 0.334;
h = 0.001;   %[m]
f = linspace(20,20000, 20000);
%% d) Propagation speed for quasi-longitudinal and longitudinal waves

c_L = sqrt(E/(rho*(1-v^2)))
c_LL = sqrt(E*(1-v)/(rho*(1+v)*(1-2*v)))

%% e) Propagation speed for bending waves
close all;
 
v_f = sqrt(1.8*f*h*c_L);

figure
semilogx(f, v_f, LineWidth=1.2)
grid minor

%% f) modal frequencies of the five bending modes of the plate

close all;

plate_facs.values = [
    0.4694, 2.08, 3.41, 5.00, 6.82 ,...
    3.89, 5.95, 8.28, 10.87, 13.71,...
    8.72, 11.75, 15.06, 18.63, 22.47
    ];

plate_facs.idx = {
    [0,1],[1,1],[2,1],[3,1],[4,1], ...
    [0,2],[1,2],[2,2],[3,2],[4,2], ...
    [0,3],[1,3],[2,3],[3,3],[4,3]};

plate_freqs = zeros(18, 1);
plate_freqs(1) = plate_facs.values(1) * c_L * h / a^2;
for i=2:length(plate_facs.values)
    plate_freqs(i) = plate_facs.values(i)*plate_freqs(1);
    %disp(['freq: ', num2str(plate_freqs(i)), ', mn: ', num2str(cell2mat(plate_facs.idx(i)))] )
end

for i=1:length(plate_facs.values)
    disp(['freq: ', num2str(plate_freqs(i)), ', mn: ', num2str(cell2mat(plate_facs.idx(i)))] )
end

%% Part 3) 
rho_vol_str = 5000; %[kg/m^3]
a_str = 0.001;
sur_str = a_str^2 * pi;
L_str = 0.4; %[m]
%% g)

string_freq = plate_freqs(1)
c = (2*string_freq)* L_str
rho_str = rho_vol_str * sur_str
T_str = c^2 * rho_str

%% h)

E_str = 200e+9;  %[Pa]
freqs_modes = zeros(5,1);
K_str = a_str/2;
B_str = (pi*E_str*sur_str*K_str^2)/(T_str*L_str^2);
for nn=1:length(freqs_modes)
    freqs_modes(nn)= nn* string_freq * sqrt(1+B_str*nn^2)*(1+(2/pi)*(B_str^(1/2))+(4/pi^2)*B_str);
end

% fn = nf ◦1√1 + Bn2[1 + 2/πB0.5 + (4/π2)B].