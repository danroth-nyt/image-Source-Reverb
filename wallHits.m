%**************************************************************************
% wallHits.m
% Adam Hill
% February 11, 2008
%
% Function to calculate how many impulses have hit a wall.
%
% Input arguments: n    =   current index number
% Ouput values: wall_1, wall_2
%
%**************************************************************************
function [wall_1 wall_2] = wallHits(n)

% Find the umber of times each wall is encountered
if n < 0
    wall_1 = ceil(abs(n)/2);                % how many wall_1 hits
    wall_2 = floor(abs(n)/2);               % how many wall_2 hits
elseif n > 0
    wall_1 = floor(n/2);                    % how many wall_1 hits
    wall_2 = ceil(n/2);                     % how many wall_2 hits
else
    wall_1 = 0;                             % direct sound
    wall_2 = 0;                             % direct sound
end

%**************************************************************************
% END OF FUNCTION
%**************************************************************************