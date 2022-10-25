%%%%%%% HOMEWORK 2 %%%%%%%%

clc; close all; clear;

addpath(genpath('HW2'))
T = 10; % N/mc 

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
for i=1:18
    subplot(3,6, i)
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

%% PART 1- c 2
Fs = 48000;
T = 1/Fs;
t_time = 2;
L = t_time*Fs; % number of samples
t = linspace(0, t_time, L);  

%f = linspace(1,200, 2000);
%len = length(f);
%t = linspace(0, 50, len);
force = 0.1*exp(-(t-0.03).^2/0.01^2);
force_omega = fft(force);
%P2 = abs(force_omega/L);
P2 = force_omega/L;
P1 = P2(1:L/2);
P1(2:end-1) = 2*P1(2:end-1);


%force_omega = force_omega(1:len/2);
%f = Fs*(0:(L/2))/L;
Fn = Fs/2; % Nyquist frequency
f = linspace(0, Fn, L/2);
% close all
figure;
plot(1000*t,force)
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')


figure
plot(f, abs(P1), LineWidth=1.2);
xlim([0,500]);

grid minor



%f = linspace(1,200, 2000);
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
        xlim([0, 200])
        end
    end
end
subplot(graphmodes, graphmodes, sub2ind([graphmodes, graphmodes], graphmodes, 1));
plot(f, zeros(1, length(f)));
grid minor
xlim([0, 200])
subplot(graphmodes, graphmodes, sub2ind([graphmodes, graphmodes], 1, graphmodes));
plot(f, zeros(1, length(f)));
grid minor
xlim([0, 200])

x = zeros(2, modes);
% phi = [15, 195];
radius = 0.075;

for nn = 1:2
    for i=1:modes
        idx = cell2mat(J.idx(i));
        m = idx(1);
        kn = freqs(i)*2*pi/c;
        %Jm = besselj(m, kn*r);
        %Jm = Jm(r==radius);
        Jm =  besselj(m, kn*radius);
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
    % zer = zeros(1, 200);
    H_tmp = diag(H(:,ii));
    H21(ii) = x(1,:) * H_tmp * x(2,:)';
    mob(ii) = 1j*2*pi*f(ii)*H21(ii);
end

figure
plot(f, abs(H21), LineWidth=1.2, LineWidth=1.2);
title("H21 2")
xlim([0,200])
grid minor



%len = 2*length(f);
%t = linspace(0, 50, len);
%force = 0.1*exp(-(t-0.03).^2/0.01^2);
%force_omega = fft(force);
%force_omega = force_omega(1:len/2);

% close all
%figure
%plot(f, abs(force_omega), LineWidth=1.2);
%grid minor

% displacement = force_omega.*mob;
displacement = P1.*H21;

figure
plot(f, abs(displacement), LineWidth=1.2, LineWidth=1.2);
title("H21")
xlim([0,300])
grid minor
% close all
%figure
%plot(f, abs(displacement), LineWidth=1.2);
%grid minor


% displacement = displacement(1:end-1);
dsp = ifft(displacement, L);
dsp = dsp*L;

% close all
figure
t2 = t(1:length(t)/2);
plot(t, real(dsp), LineWidth=1.2);
hold on
plot(t, force*3e-5, LineWidth=1.2);
grid minor
xlim([0, 1])




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

%%
sound(real(dsp)/max(real(dsp)),Fs)



%% Part 2) 
E = 69e+9;  %[Pa]
rho = 2700; %[kg/m^3]
v= 0.334;
h = 0.001;   %[m]
f = linspace(20,20000, 20000);
M_plate = rho*h*pi*a^2;
%% d) Propagation speed for quasi-longitudinal and longitudinal waves

c_L = sqrt(E/(rho*(1-v^2)))
c_LL = sqrt(E*(1-v)/(rho*(1+v)*(1-2*v)))

%% e) Propagation speed for bending waves
close all;
 
v_f = sqrt(1.8*f*h*c_L);

figure
semilogx(f, v_f, LineWidth=1.2)
grid minor
title('Bending waves propagation speed as function of the frequency', Interpreter='tex');
xlabel('Frequency [Hz]',Interpreter='tex')
ylabel('v(f) [m/s]',Interpreter='tex')

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

plate_freqs = zeros(15, 1);
plate_freqs(1,1) = plate_facs.values(1) * c_L * h / a^2;

for i=2:length(plate_facs.values)
    plate_freqs(i) = plate_facs.values(i)*plate_freqs(1);
    %plate_freqs(i,2) = i;
    %disp(['freq: ', num2str(plate_freqs(i)), ', mn: ', num2str(cell2mat(plate_facs.idx(i)))] )
end
[plate_freqs,plate_order] = sort(plate_freqs);
for i=1:length(plate_facs.values)
    disp(['freq: ', num2str(plate_freqs(i)), ', mn: ', num2str(cell2mat(plate_facs.idx(plate_order(i))))] )
end
plate_freqs2 = plate_freqs(1:5);
%% Part 3) 
rho_vol_str = 5000; %[kg/m^3]
a_str = 0.001;
sur_str = a_str^2 * pi;
L_str = 0.4; %[m]
%% g)

string_freq = plate_freqs(1)
c_str = (2*string_freq)* L_str
mu_str = rho_vol_str * sur_str
T_str = c_str^2 * mu_str

%% h) frequencies of the first five modes of the string with clamped edges 

E_str = 200e+9;  %[Pa]
freqs_modes_clamped = zeros(5,1);
K_str = a_str/2;
B_str = (pi*E_str*sur_str*K_str^2)/(T_str*L_str^2);
for nn=1:length(freqs_modes_clamped)
    freqs_modes_clamped(nn)= nn* string_freq * sqrt(1+B_str*nn^2)*(1+(2/pi)*(B_str^(1/2))+(4/pi^2)*B_str);
end

% fn = nf ◦1√1 + Bn2[1 + 2/πB0.5 + (4/π2)B].

%% h) frequencies of the first five modes of the string with supported edges 
E_str = 200e+9;  %[Pa]
freqs_modes_supported = zeros(5,1);
K_str = a_str/2;
B_str = (pi*E_str*sur_str*K_str^2)/(T_str*L_str^2);
for nn=1:length(freqs_modes_supported)
    freqs_modes_supported(nn)= nn* string_freq * sqrt(1+B_str*nn^2);
end


%% i) 
Q_b = 50;
numOfHarms = 5;
string_freqs_noDamp = zeros(numOfHarms,1);
for nn = 1:numOfHarms 
   string_freqs_noDamp(nn) = (nn*c_str*pi)/(2*pi*L_str);
end 

m_str = L_str* mu_str;
plate_string_freqs = zeros(5,4);
plate_string_freqs(:,1) = plate_freqs2;
plate_string_freqs(:,2) = freqs_modes_supported;
plate_string_freqs(:,3) = freqs_modes_clamped;
plate_string_freqs(:,4) = string_freqs_noDamp;
plate_freqs2 = round(plate_freqs2,4)
string_freqs_noDamp = round(string_freqs_noDamp,4)

% first column: plate resonance frequencies
% second column: string not stiff resonance frequencies
% third column: 0 weak coupling, 1 strong coupling with same frequency case

% fourth column: m/n^2*M * 10^3
% fifth column: delta frequencies;
% sixth column: first splitted freq, omega positive;
% seventh column: second splitted freq, omega negative;
% eighth column: ws-wb/wb
plate_string_coupling=zeros(5,8); 
plate_string_coupling(:,1) = plate_freqs2;
plate_string_coupling(:,2) = string_freqs_noDamp;
for i= 1:length(string_freqs_noDamp)
    if(string_freqs_noDamp(i)==plate_freqs2(i))
        if(m_str*(4*Q_b^2)/(i^2*M_plate)>pi^2)
            plate_string_coupling(i,3) = 1;   
        end
        
    end
    plate_string_coupling(i,4) = m_str/(i^2*M_plate)*10^3;
end    

% done by hand cause I don't have a fucking matlab curve
plate_string_coupling(1,5) = 0.078;
plate_string_coupling(5,5) = 0;
for i= 1:length(string_freqs_noDamp)
    if plate_string_coupling(i,5) ~= 0
        plate_string_coupling(i,6) = (1+plate_string_coupling(i,5))...
            *plate_string_coupling(i,2);
        plate_string_coupling(i,7) = (1-plate_string_coupling(i,5))...
            *plate_string_coupling(i,2);
        
    else
        if plate_string_coupling(i,3) ~= 1
            plate_string_coupling(i,8) = (plate_string_coupling(i,2)-plate_string_coupling(i,1))/plate_string_coupling(i,1);
        end

    end

end   

<<<<<<< HEAD
plate_string_coupling(2,9) = 0.01;
plate_string_coupling(2,10) = -0.05;
plate_string_coupling(3,9) = 0;
plate_string_coupling(3,10) = -0.1202;
plate_string_coupling(4,9) = 0.04;
plate_string_coupling(4,10) = -0.01;

plate_string_coupling(5,6) = plate_string_coupling(5,1);
plate_string_coupling(5,7) = plate_string_coupling(5,1);


for i=1:length(string_freqs_noDamp)
    if plate_string_coupling(i,8) ~=0
        plate_string_coupling(i,6) = (1+plate_string_coupling(i,9))*plate_string_coupling(i,1);
        plate_string_coupling(i,7) = (1+plate_string_coupling(i,10))*plate_string_coupling(i,1);
    end
end
%plate_string_coupling_omega = plate_string_coupling*2*pi
=======
>>>>>>> parent of 4e936a0 (added coupled system freqs calculations )

% m/(n2M) < π2/(4Q2B),
