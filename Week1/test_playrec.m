% Tests the functionality of initparams() and recplay.mdl

%% Generate a sinewave
clear all
fs = 16000; % Sampling rate [Hz]
f0 = 1500; % Frequency of sinewave [Hz]
t_max = 2 ; % Length of the signal [s]


% Generate sinewave
duration = 2;
t = linspace(0,duration,duration*fs); %duration * fs = aantal samples

sinewave = sin(2*pi*f0*t)';

% Call to initparams()
%function [ simin,nbsecs,fs ] = initparams( toplay, fs, varargin )
[simin,nbsecs,fs] = initparams(sinewave,fs,3);

%% Play and record signal
% Call to recplay.mdl to play simin and record simout
sim('recplay'); 

% Retrieve recorded output
out=simout.signals.values; 

%% Playback
% Play the sound using sound() or soundsc()
sound(out)