%%%%%%% HOMEWORK 2 %%%%%%%%

clc; close all; clear;

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

f = linspace(1,200, 2000);
omega = 2*pi*f;
Q = 25;
modes = 18;

close all
figure
H = zeros(18, length(f));

for i = 1:modes
    omega0i = (2*pi*freqs(i));
    H(i, :) = 1./(-omega.^2 + 1j*omega*omega0i/Q + omega0i^2);
    
    subplot(modes,modes, sub2ind([modes,modes], i,i));
    plot(f, abs(H(i,:)));
end

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
plot(f, abs(mob), LineWidth=1.2)
grid minor

len = 2*length(f);
t = linspace(0, 50, len);
force = 0.1*exp(-(t-0.03).^2/0.01^2);
force_omega = fft(force);
force_omega = force_omega(1:len/2);

% close all
figure
plot(f, abs(force_omega));

% displacement = force_omega.*mob;
displacement = force_omega.*H21;

% close all
figure
plot(f, abs(displacement));


dsp = ifft(displacement, len);

% close all
figure
plot(t, real(dsp))
xlim([0, 9])
grid minor



%% Part 2) 
E = 69e+9 %[Pa]
rho = 2700 %[kg/m^3]
v= 0.334
%% d)

c_L = 5


