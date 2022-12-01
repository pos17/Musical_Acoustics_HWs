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
clear all;
close all;
clc;

addpath('Functions')
simulink_folder = './';       % Simulink projects folder
addpath(simulink_folder);

%% Setup
fs = 44100;                         % Sampling frequency
signalLen = 5;                      % Signal length
t = [0:1/fs:signalLen-1/fs];        % Time axis

fileName = 'GuitarModel.wav';     % Audio file path

%% Simulation
% run the simulink simulation using the command sim (see doc sim).
sim('GuitarModel2'); 

% The variable I contains non constant time intervals between samples.
% Resample the data using resample function in order to obtain an equally
% sampled signal I1
I1 = resample(I, t);
%%
% Plot the resampled signal in time
figure(1)
plot(I1.Time, I1.Data);

% Normalize the signal
soundWave = I1.data;
soundWave = soundWave./max(abs(soundWave));

%% Plot and play
% Plot the signal frequency content as magnitude and phase
[S, magS, angleS, f, df] = myFFT(soundWave, fs);

h = figure(2);
plotFFT_linearFreqScale(magS, angleS, f, df, fs, 2000, h);

figure(3);
plot(t, soundWave);
title('Signal in time'), xlabel('Time [s]')

sound(soundWave, fs);                           % Play the sound

disp('Save file on disk...')                    % Save on disk
audiowrite(fileName, soundWave, fs);

%% Body and string 

fileName = 'GuitarModelString.wav';

%% Simulation
% run the simulink simulation using the command sim (see doc sim).
sim('GuitarModel3'); 

% The variable I contains non constant time intervals between samples.
% Resample the data using resample function in order to obtain an equally
% sampled signal I1
I2 = resample(I_2, t);
%%
% Plot the resampled signal in time
figure(4)
plot(I2.Time, I2.Data);

% Normalize the signal 

soundWave2 = I2.data;
soundWave2 = soundWave2./max(abs(soundWave2));

%% Plot and play
% Plot the signal frequency content as magnitude and phase
[S, magS, angleS, f, df] = myFFT(soundWave2, fs);

h = figure(5);
plotFFT_linearFreqScale(magS, angleS, f, df, fs, 2000, h);

figure(6);
plot(t, soundWave2);
title('Signal in time'), xlabel('Time [s]')

sound(soundWave2, fs);                           % Play the sound

disp('Save file on disk...')                    % Save on disk
audiowrite(fileName, soundWave2, fs);
