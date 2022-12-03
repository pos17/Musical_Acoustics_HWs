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
simulink_folder = './PreviousVersion/';       % Simulink projects folder
addpath(simulink_folder);

%% Setup
fs = 48000;                         % Sampling frequency
signalLen = 5;                      % Signal length
% t = [0:1/fs:signalLen-1/fs];        % Time axis
t = linspace(0, signalLen, (signalLen*fs)+1);

fileName = 'GuitarModel.wav';     % Audio file path


%% Simulation
% run the simulink simulation using the command sim (see doc sim).
flag = 0;
sim('GuitarModel2'); 
%%

% The variable I contains non constant time intervals between samples.
% Resample the data using resample function in order to obtain an equally
% sampled signal I1
% I1 = resample(I, t);

figure();
plot(I.Time, I.Data);
close all
H = fft(I.Data)./fft(V.Data);
H = H(1:floor(length(H)/2));

frq = linspace(0,fs/2, length(H));

% g = figure();
g = figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
FRF_Plot('H', H, frq, g);
sgtitle("Admittance of the Guitar Body ")


[pks, locs] = findpeaks(abs(H));

resonances = frq(locs)';
clc
disp("resonances of the system");
disp(resonances);
% %%
% close all
% figure(15)
% bode(H)
% 
% V1 = resample(V, t);
%%
% Plot the resampled signal in time
flag = 2;
sim('GuitarModel2'); 
%%
close all
figure('Renderer', 'painters', 'Position', [100 100 1000 600])
plot(I.Time, I.Data);
grid minor
xlabel('time [s]'); ylabel('current [A]');
title('Velocity of top plate');


figure('Renderer', 'painters', 'Position', [100 100 1000 600])
plot(V.Time, V.Data);

% Normalize the signal
soundWave = I.data;
soundWave = soundWave./max(abs(soundWave));

%% Plot and play
% Plot the signal frequency content as magnitude and phase
% [S, magS, angleS, f, df] = myFFT(soundWave, fs);
% [V_f, magV, angleV, f, df] = myFFT(V1.data, fs);

close all
h = figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
S = fft(I.Data);
S = S(1:floor(length(S)/2));

FRF_Plot(' . ', H, frq, h);
hold on
FRF_Plot('' ,S, frq, h);
legend("admittance", "top plate vel")
sgtitle("Top Plate Velocity vs Admittance")

% j = figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
% F = fft(V.Data);
% F = F(1:floor(length(F)/2));
% 
% FRF_Plot('Force' ,F, frq, j);
% sgtitle("Frequency response of the exciter force")


% plotFFT_linearFreqScale(magS, angleS, f, df, fs, fs/4, h);
% g = figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
% plotFFT_linearFreqScale(magV, angleV, f, df, fs, fs/4, g);

% j = figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
% plotFFT_linearFreqScale(magS./magV, angleS-angleV, f, df, fs, fs/4, j);



%%
figure(3);
plot(t, soundWave);
title('Signal in time'), xlabel('Time [s]')

sound(soundWave, fs);                           % Play the sound

disp('Save file on disk...')                    % Save on disk
audiowrite(fileName, soundWave, fs);


% H = abs(S./V_f);
% g = figure();
% plotFFT_linearFreqScale(H, angle(S./V_f), f, df, fs, 2000, g);


%% Body and string 

fileName = 'GuitarModelString.wav';

%% Simulation
% run the simulink simulation using the command sim (see doc sim).
sim('GuitarModel3'); 

% The variable I contains non constant time intervals between samples.
% Resample the data using resample function in order to obtain an equally
% sampled signal I1

I2 = resample(I_2, t);
% I2 = I_2;

% Plot the resampled signal in time
% figure(4)
% plot(I2);
% figure(12)
% spectrogram(I_2.Data, hamming(128), 64);

% Normalize the signal 

soundWave2 = I2.data;
soundWave2 = soundWave2./max(abs(soundWave2));

%% Plot and play
% Plot the signal frequency content as magnitude and phase
S = fft(I2.Data); 
S = S(1:floor(length(S)/2));
jj = figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
frq = linspace(0,fs/2, length(S));
FRF_Plot('S', S, frq, jj)

figure();
plot(t, I2.Data);
title('Top Plate Velocity Signal');
xlabel('Time [s]'); ylabel('Current [A]');
grid minor
%%
close all
figure()
spectrogram(soundWave2, blackman(2048), 64, 1024, fs, 'psd');
xlim([0, 5]);
%% PRINTING FILE

sound(soundWave2, fs);                          % Play the sound
%%
disp('Save file on disk...')                    % Save on disk
audiowrite(fileName, soundWave2, fs);
