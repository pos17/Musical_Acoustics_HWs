%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling of musical instruments homework.                               %
% Numerical simulation of a guitar                                        %
% Physical model for a guitar exploiting electric analogs.                %
%                                                                         %
% Musical Acoustics course                                                %
% Ostan Paolo                                                             %
% Donà Stefano                                                            %
% 2022                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
% in this file we didn't use the functions myFFT() and 
% plotFFT_linearFreqScale() because we noticed that the fft returned by the
% first function wasn't smooth enough in order to collect qualitative and
% quantitave informations from it. Instead, the built-in fft of MATLAB gave
% us the smoothness we needed. Moreover, given the high differences in
% terms of order of magnintude of the different peaks we built our own
% function FRF_Plot that is an easier modification of the second given
% function, and it plots the magnitude in dB scale in order to understand
% better the behavior of the system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;
close all;
clc;

addpath('Functions')
simulink_folder = './PreviousVersion/';       % Simulink projects folder
addpath(simulink_folder);

if not(isfolder("plots\"))
    mkdir("plots\")
end


%% Setup
fs = 48000;                         % Sampling frequency
signalLen = 5;                      % Signal length    
t = linspace(0, signalLen, (signalLen*fs));  % Time axis


%% SIMULATION OF IMPULSE RESPONSE
% run the simulink simulation using the command sim (see doc sim).
flag = 0; % commutation of switch to the impulse excitation
sim('GuitarModel2'); 


%% PLOTTING OF THE ADMITTANCE OF THE GUITAR BODY

% The variable I contains non constant time intervals between samples.
% Resample the data using resample function in order to obtain an equally
% sampled signal I1
I = resample(I_out, t);
V = resample(V_out, t);

H = fft(I.Data)./fft(V.Data);
H = H(1:floor(length(H)/2));

frq = linspace(0,fs/2, length(H));

g = figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
FRF_Plot('H', H, frq, 2000, g);
sgtitle("Admittance of the Guitar Body ")
% saving plot
saveas(gcf, ".\plots\GTR_Ex1_Admittance.png");


[pks, locs] = findpeaks(abs(H));
resonances = frq(locs)';

disp("resonances of the system");
disp(resonances);


%% SIMULATION OF HARMONIC EXCITATION
flag = 2; % commutation of the switch to the damped squared signal
sim('GuitarModel2'); 
I = resample(I_out, t);


%% PLOTTING THE 

figure('Renderer', 'painters', 'Position', [100 100 1000 600])
plot(I.Time, I.Data);
grid minor
xlabel('time [s]'); ylabel('current [A]');
title('Velocity of top plate');
% saving plot
saveas(gcf, ".\plots\GTR_Ex1_VelTime.png");

figure('Renderer', 'painters', 'Position', [100 100 1000 600])
plot(I.Time, I.Data);
grid minor
xlabel('time [s]'); ylabel('current [A]');
title('Velocity of top plate zoomed');
xlim([0 0.1])
% saving plot
saveas(gcf, ".\plots\GTR_Ex1_VelTimeZoomed.png");


%% Plot and play
% Plot the signal frequency content as magnitude and phase
S = fft(I.Data);
S = S(1:floor(length(S)/2));

h = figure('Renderer', 'painters', 'Position', [100 100 1000 600]);

FRF_Plot('', H, frq, 2000, h);
hold on
FRF_Plot(' . ' , S, frq, 2000, h);
legend(["admittance", "top plate vel"], Position=[0.790099999227524,0.47358333269755,0.115800001544952,0.060833334604899])
sgtitle("Top Plate Velocity vs Admittance")
% saving plot
saveas(gcf, ".\plots\GTR_Ex1_AdmVSVel.png");


%% LISTENING AND SAVING SOUNDWAVE

fileName = 'GuitarModelNoString.wav';     % Audio file path

% Normalize the signal
soundWave = I.Data;
soundWave = soundWave./max(abs(soundWave));

sound(soundWave, fs);                           % Play the sound

disp('Save file on disk...')                    % Save on disk
audiowrite(fileName, soundWave, fs);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   SIMULATION OF BODY AND REAL STRING                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Body and string Simulation
% run the simulink simulation using the command sim (see doc sim).
sim('GuitarModel3'); 

% The variable I contains non constant time intervals between samples.
% Resample the data using resample function in order to obtain an equally
% sampled signal I_2

I2 = resample(I_2, t);

% Plot the resampled signal in time
figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
plot(I2);
xlabel('time [s]'); ylabel('current [A]')
title('Velocity of top plate w/ string model')
grid minor
% saving plot
saveas(gcf, ".\plots\GTR_Ex2_VelTime.png");

% Normalize the signal 
soundWave2 = I2.Data;
soundWave2 = soundWave2./max(abs(soundWave2));


%% PLOTTING FFT

% Plot the signal frequency content as magnitude and phase
S = fft(I2.Data); 
S = S(1:floor(length(S)/2));
jj = figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
frq = linspace(0,fs/2, length(S));
FRF_Plot('Vel', S, frq, 5500, jj)
sgtitle('Frequency response of the velocity w/ string model')
allAxis=get(gcf, 'Children');

f0 = 329.5;
xline(f0, 'k--', 'LineWidth', 2,'Parent', allAxis(3));
xline([f0*5 f0*10 f0*15], 'r--', 'LineWidth', 2,'Parent', allAxis(3));

% saving plot
saveas(gcf, ".\plots\GTR_Ex2_VelFreq.png");


%% LISTENING AND PRINTING FILE

fileName = 'GuitarModelString.wav';

sound(soundWave2, fs);                          % Play the sound

disp('Save file on disk...')                    % Save on disk
audiowrite(fileName, soundWave2, fs);
