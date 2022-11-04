
close all;
clear all;
clc;

% Simulation Parameters
Fs = 2000;                  % Sampling Frequency
signalLen = 3;              % Simulation Duration
z0 = 2.224247e+05;  % Filter characteristic impedance


%% ------------ Filter roll-off increasing number of elements ---------- %%

set_param('exercise_4a_complete/R1', 'R', num2str(z0));
set_param('exercise_4a_complete/R2', 'R', num2str(z0));

figure; hold on;

label = ['a'; 'b'; 'c'];
for i = 1:3

    % Load Simscape model
    disp(['run model ' label(i)]);
    open_system(['exercise_4' label(i) '_complete.slx'], 'loadonly');
    
    % Perform Simulation
    if strcmp(label(i), 'a')
            set_param('exercise_4a_complete/R1', 'R', num2str(z0));
    end
    simulation = sim(['exercise_4' label(i) '_complete.slx'], signalLen);
    
    % Compute Frequency Response
    input = simulation.simin.data;
    output = simulation.simout.data;
    f = 0:Fs/length(input):Fs-(1/length(input));    
    H = db(abs(fft(output) ./ fft(input)));
    
    % Plot
    plot(f, H);
    pause(0.05);    

end

hold off;
title('Low-Pass Muffler - Frequency Response');
xlabel('Frequency [Hz]');
xlim([0 500]);
ylabel('Magnitude [dB]');
legend('1 section', '2 sections', '3 sections');


%% ---------------- Filter Response varying source load ---------------- %%

open_system('exercise_4a_complete.slx', 'loadonly');
set_param('exercise_4a_complete/R1', 'R', num2str(z0));
set_param('exercise_4a_complete/R2', 'R', num2str(z0));

figure; hold on;
for i = 0:4
    
    % Change Source Load
    temp_R = z0/(2^i);
    set_param('exercise_4a_complete/R1', 'R', num2str(temp_R));
    
    % Simulation
    disp(['Run model with R1 = ' num2str(temp_R)]);
    simulation = sim(['exercise_4a_complete.slx'], signalLen);
    
    % Compute Frequency Response
    input = simulation.simin.data;
    output = simulation.simout.data;
    f = 0:Fs/length(input):Fs-(1/length(input));    
    H = db(abs(fft(output) ./ fft(input)));
    
    % Plot
    plot(f, H);
    pause(0.05);

end

hold off;
title('Low-Pass Muffler - Frequency Response');
xlabel('Frequency [Hz]');
xlim([0 500]);
ylabel('Magnitude [dB]');
legend('z0', 'z_0/2', 'z_0/4', 'z_0/8', 'z_0/16');


%% ---------------- Filter Response varying output load ---------------- %%

set_param('exercise_4a_complete/R1', 'R', num2str(z0));
set_param('exercise_4a_complete/R2', 'R', num2str(z0));

figure; hold on;
for i = 0:4
    
    % Change Source Load
    temp_R = z0/(2^i);
    set_param('exercise_4a_complete/R2', 'R', num2str(temp_R));
    
    % Simulation
    disp(['Run model with R1 = ' num2str(temp_R)]);
    simulation = sim(['exercise_4a_complete.slx'], signalLen);
    
    % Compute Frequency Response
    input = simulation.simin.data;
    output = simulation.simout.data;
    f = 0:Fs/length(input):Fs-(1/length(input));    
    H = db(abs(fft(output) ./ fft(input)));
    
    % Plot
    plot(f, H);
    pause(0.05);

end

hold off;
title('Low-Pass Muffler - Frequency Response');
xlabel('Frequency [Hz]');
xlim([0 500]);
ylabel('Magnitude [dB]');
legend('z0', 'z_0/2', 'z_0/4', 'z_0/8', 'z_0/16');