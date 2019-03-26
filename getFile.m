% getFile.m
% Adam Hill
% June 24, 2008
%
% Function to get a user defined .wav file for use in responseGUI.m
%
% Output: .wav file in vector form

function [y, SR] = getFile(label)

% get the desired file choice from the user
[filename, pathname] = uigetfile('*.wav', strcat('Pick a .wav file for: ', label));

% Check is a file was selected, convert .wav to vector if it was
if isequal(filename,0)
   y = [];
else
   [y, SR] = audioread(fullfile(pathname, filename));
   y = y(:, 1)';
end