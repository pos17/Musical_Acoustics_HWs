%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling of musical instruments homework.                               %
% Numerical simulation of piano strings                                   %
% Physical model for a struck string using finite difference.             %
%                                                                         %
% Musical Acoustics course                                                %
% Mirco Pezzoli                                                           %
% 2022                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%% Setup
% define the string parameters and the simulation variables defined
% according to the provided values and the numerical implementation.
% We want to implement the finite difference scheme of a piano string tuned
% as C2.

% Temporal sampling parameters
Fs = 16000;                     % [Hz] Following the paper
%Fs = 4*44100;                  % [Hz] Following the slides
%Fs = 44100;                  % [Hz] Following the slides

T = 1/Fs;
dur = 8;                        % [s]
N = dur*Fs;
%t = linspace(0, Fs, N);


% Boundary          

z_l = 1e+20;
z_b = 1000;

% String parameters
f_1 = 65.4;                     % [Hz]  Fundamental note
str_L = 1.92;                   % [m]   string length
M_s = 35e-3;                    % [Kg]  String mass
%T_e = 750;                     % [N]   string tension given  by the slides
 

b_1 = 0.5;                      %       air damping coefficient
b_2 = 6.25e-9;                  %       string internal friction coefficient
epsilon = 7.5 * 10^(-6);        %       string stiffness parameter
k = epsilon;                    %       string stiffness coefficient

rho = M_s / str_L;              % [Kg/m]string linear density    

T_e = 4*str_L^2*rho*f_1^2;      % [N]   string tension calculated  
c = sqrt(T_e/rho);              % [m/s] string propagation velocity

% Spatial sampling parameters
%M = X = T * c       % maximum value which allows CFL conditions
%M = 521;            % fixed by the slides



% Aliasing condition


% Number of maximum spatial steps
gamma = Fs/(2*f_1);
M_max = sqrt((-1+sqrt(1+16*epsilon*gamma^2))/(8*epsilon));
X_max = sqrt((1/2)*(c^2*T^2 + 4*b_2*T + sqrt(((c^2)*T^2 + 4*b_2*T)^2 + 16*k^2*T^2)));
X_maxCFL = c/Fs;
%M_max = str_L/X_max;
M = floor(M_max);
%M=120;
X = str_L/M;

% Integer values
lambda = c*T/X;
mu = k^2/(c^2*X^2);
v= (2*b_2*T)/(X^2);


% Spatial sampling



% FD parameters
a_1 = (-lambda^2*mu)/(1+b_1*T);
a_2 = (lambda^2 + (4*lambda^2*mu) + v)/(1+b_1*T);
a_3 = (2-2*lambda^2-6*lambda^2*mu-2*v)/(1+b_1*T);
a_4 = (-1+b_1*T+2*v)/(1+b_1*T);
a_5 = (-v)/(1+b_1*T);
a_F = (T^2/rho)/(1+b_1*T);

% Hammer parameters
M_H = 4.9e-3;                   % [Kg]  mass of the hammer
w = 0.2;                        %       width of the hammer spatial window 𝑔
p = 2.3;
Vh_0 = 2.5;                     % [m/s] initial hammer velocity
b_H = 1e-4;                     % [1/s] fluid damping coefficient
K = 4e+8;                       %       hammer felt stiffness
a=0.12;                         %       x_0/L relative position of hammer strike
a_M = floor(M*a);               %       relative hammer position in space samples

% Hammer contact window definition
w_M= round(w/X);
g = hann(w_M);

%PDE Coefficients:

% Bridge boundary coefficients
b_R1 = (2-2*lambda^2*mu-2*lambda^2)/(1+b_1*T+z_b*lambda);
b_R2 = (4*lambda^2*mu+2*lambda^2)/(1+b_1*T+z_b*lambda);
b_R3 = (-2*lambda^2*mu)/(1+b_1*T+z_b*lambda);
b_R4 = (-1+b_1*T+z_b*lambda)/(1+b_1*T+z_b*lambda);
b_RF = (T^2/rho)/(1+b_1*T+z_b*lambda);

% Left hand (hinged string end) boundary coefficients
b_L1 = (2-2*lambda^2*mu-2*lambda^2)/(1+b_1*T+z_l*lambda);
b_L2 = (4*lambda^2*mu+2*lambda^2)/(1+b_1*T+z_l*lambda);
b_L3 = (-2*lambda^2*mu)/(1+b_1*T+z_l*lambda);
b_L4 = (-1+b_1*T+z_b*lambda)/(1+b_1*T+z_l*lambda);
b_LF = (T^2/rho)/(1+b_1*T+z_l*lambda);

% Hammer felt parameters

d_1 = 2/(1+(b_H*T)/(2*M_H));
d_2 =(-1 + (b_H*T)/(2*M_H))/(1 + (b_H*T)/(2*M_H));
d_F = (-T^2/M_H)/(1+(b_H*T)/(2*M_H));


%% Computation of the FD scheme
% Initialization

% y string displacement vector 
y=zeros([N,M]);     % each row is a time instant considered. 
                    % each column is a space fragment

av_y = zeros(N,1);  % averaged displacement over 12 spatial samples

% F_H time vector
F_H=zeros(N,1);     % each row is a time instant considered. 
                    % each column is a space fragment

F=zeros(N,M);      % 

space = linspace(0, str_L, M);
j_start = a_M-(round(w_M/2));
disp(j_start)
for i=1:N
    for j= 1:size(g)
        F(i,(j)+j_start) = g(j); % dovrebbe iniziare da 0 l'offset??? BOH

    end

end


%% g(x,x0) force density
figure('Renderer', 'painters', 'Position', [100 100 1000 500])
%plot(space/str_L, F(10,:), LineWidth=1.2);
plot(space, F(10,:), LineWidth=1.2);

grid minor;
% surf(F);
disp(F(10,a_M));
xlabel("x [m]")
ylabel("g ")
title('$g(x,x_{0})$ shape and position of the force density function','interpreter','latex','FontSize',15)

delete("./pianoAssets/ForceDensityShape.png");
saveas(gcf, "./pianoAssets/ForceDensityShape.png");
%%  LOOPS CALCULATIONS 

% eta hammer displacement vector 
eta = zeros(N,1);
eta(1) = 0;         % eta time= 0
eta(2) = Vh_0*T; 

F_H(1) = K* abs(eta(1)-y(1,a_M))^p;
F_H(2) = K* abs(eta(2)-y(2,a_M))^p;
% Computation loop



for in=2:N-1          %   time loop

    % force time computation
   
    for im=1:M      %   space loop

        % solving displacement boundary conditions
        switch im
            % m=0
            case 1
                y(in+1,im) = b_L1*y(in,im)+b_L2*y(in,im+1)+b_L3*y(in,im+2)...
                    +b_L4*y(in-1,im)+b_LF*F(in,im);
            % m=1
            case 2
                y(in+1,im) = a_1*(y(in,im+2)-y(in,im)+2*y(in,im-1))+a_2*(y(in,im+1)+y(in,im-1))+ ...
                    +a_3*(y(in,im))+a_4*(y(in-1,im))+a_5*(y(in-1,im+1)+y(in-1,im-1))+a_F*F(in,im);
            % m=M-1
            case M-1
                y(in+1,im) = a_1*(2*y(in,im+1)-y(in,im)+y(in,im-2))+a_2*(y(in,im+1)+y(in,im-1))+ ...
                    +a_3*(y(in,im))+a_4*(y(in-1,im))+a_5*(y(in-1,im+1)+y(in-1,im-1))+a_F*F(in,im);
            % m=M
            case M
                y(in+1,im)= b_R1*y(in,im)+b_R2*y(in,im-1)+b_R3*y(in,im-2)+b_R4*y(in-1,im)+b_RF*F(in,im);
            otherwise
                y(in+1,im)= a_1*(y(in,im+2)+y(in,im-2))+a_2*(y(in,im+1)+y(in,im-1))+a_3*y(in,im)+a_4*(y(in-1,im))+a_5*(y(in-1,im+1)+y(in-1,im-1))+a_F*(F(in,im));
        end
    end

    eta(in+1)= d_1*eta(in)+d_2*eta(in-1)+d_F*F_H(in);

    if eta(in+1)>= y(in+1,a_M)
            F_H(in+1) = K* abs(eta(in+1)-y(in+1,a_M))^p;
            F(in+1,:)=F(in+1,:)*F_H(in+1);
    else
            F_H(in+1) =0;
            F(in+1,:) = 0;
    end

    
    for im=-5:6                     % NOT POSSIBLE TO AVERAGE IT PERFECTLY
                                    % OVER 12 SAMPLES
        av_y(in,1) = av_y(in,1)+y(in,im+a_M);
    end
    av_y(in,1) = av_y(in,1)/12;
    
    
end

%% Plot the displacement of the whole string in time
figure('Renderer', 'painters', 'Position', [100 100 1000 500])
for i = 1:(20*dur)
    iTime = round(i*1/20*Fs);

    %if rem(in,200)==1
    
    grid minor;
    plot(y(iTime,:),LineWidth=1.2);
    ylim([-0.0005 0.0005])
    grid on;
    title(strcat('t = ',num2str(i*0.05)));

    drawnow;
    pause(1/20);
    %end
end
%% Saving some samples of the string movement for report representation 
upper = 8;
lower = 0.001;
numOfPics = 9;
picTime= logspace(1,3,numOfPics);
picTime_01 = picTime/max(picTime) - min(picTime)/max(picTime)
% re-scale to be between your limits (i.e. 1 and 1.05)
picTime_08 = picTime_01*(upper-lower) + lower


for i = 1:(numOfPics)
    figure('Renderer', 'painters', 'Position', [100 100 800 500])
    picTime_08Samp = round(picTime_08*Fs);

    %if rem(in,200)==1
    
    grid minor;
    plot(y(picTime_08Samp(i),:),LineWidth=1.8);
    ylim([-0.0005 0.0005])
    grid on;
    title(strcat('string Vertical displacement y at time $t = ',num2str(picTime_08(i)),"$,",num2str(M) ," spatial samples"),'interpreter','latex','FontSize',25) 
    % title(strcat('t = ',num2str(picTime_08(i))));
    picName = strcat("./pianoAssets/PianoStringDispTime_",num2str(M),"_",num2str(i),".png");
    delete(picName);
    saveas(gcf, picName);
end

%% FFT of the averaged signal
f=linspace(0,Fs,N);
freqs = abs(fft(av_y(:,1)));
figure('Renderer', 'painters', 'Position', [100 100 1000 500])
grid minor;
plot(f,freqs, LineWidth=1.2);
xlim([0,Fs/2])


%% Plot the displacement in time
t=linspace(0, dur, N);
figure('Renderer', 'painters', 'Position', [100 100 1000 500])
plot(t, av_y(:,1), LineWidth=1.2);
grid minor

xlabel("x [m]")
ylabel("y [m]")
title('Vertical displacement of the mean of 12 samples of the string around center hammering position','interpreter','latex','FontSize',15)

delete("./pianoAssets/MeanDisplacement.png");
saveas(gcf, "./pianoAssets/MeanDisplacement.png");

%% Plot the force in time 
t=linspace(0, dur, N);
figure('Renderer', 'painters', 'Position', [100 100 1000 500])
plot(t, F(:,a_M), LineWidth=1.2);
xlim([0,0.01]);

%% Plot the force in time 
t=linspace(0, dur, N);
figure('Renderer', 'painters', 'Position', [100 100 1000 500])
plot(t, F_H(:,1), LineWidth=1.2);
xlim([0,0.1]);

%% Plot the synthesized signal play it and save it on the disk

% Play the sound
av_y_save = av_y./max(abs(av_y));
sound(av_y_save,Fs)
% pause(10);
% sound(diff(av_y_save)*10,Fs)

% Save on disk

filenamePaolo = '10868276_Ostan_piano.wav';
filenameStefano = "10868267_Donà_piano.wav";
audiowrite(filenamePaolo,av_y_save,Fs);
audiowrite(filenameStefano,av_y_save,Fs);
%clear y_av Fs


%% Spectrogram of the audio signal 

figure('Renderer', 'painters', 'Position', [100 100 800 600]);
spectrogram(av_y_save,blackman(4096),64, 4096,Fs,"power");
xlim([0,4])
title('Spectrogram of the mean string displacement around hammer center position ','interpreter','latex','FontSize',15);

% plot saving 
delete("./pianoAssets/MeanDisplacementSpectrogram.png");
saveas(gcf, "./pianoAssets/MeanDisplacementSpectrogram.png");






