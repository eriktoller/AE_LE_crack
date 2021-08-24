function [] = Plot_line(zfrom,zto,plot_prop)
% PLOT_LINE Plotting a line between two complex points points
%
%   VARIABLES
%   zfrom - start point of line (complex)
%   zto - end point of line (complex)
%   plot_prop - plotting porperties (string) given as e.g. 'black' or '.:'
%
%   LATEST UPDATE
%   2021-05-18

% Get the x and y vectors
x = [real(zfrom), real(zto)];
y = [imag(zfrom), imag(zto)];

% Plot the line
plot(x,y,plot_prop)
end

