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

f_1 = 65.4 %[Hz] % Fundamental note
% Boundary          

% String parameters

b_1 = 0.5; %air damping coefficient
b_2 = 0.5; %string internal friction coefficient
w = 0.2; % width of the hammer spatial window 𝑔
Vh_0 = 2.5 % [m/s]  initial hammer velocity
epsilon = 7.5 * 10^(-6); %string stiffness parameter
k = epsilon; % string stiffness coefficient



% Spatial sampling parameters
% Aliasing condition
% Number of maximum spatial steps

% Integer values
% Spatial sampling



% FD parameters

% Hammer parameters

% Hammer contact window definition


%PDE Coefficients:

% Bridge boundary coefficients

% Left hand (hinged string end) boundary coefficients

% Hammer felt parameters

%% Computation of the FD scheme
% Initialization

% Computation loop
%% Plot the displacement in time

%% Plot the synthesized signal play it and save it on the disk

% Play the sound



% Save on disk










