function [] = Plot_line(zfrom,zto,plot_prop)
% Plotting a line between two points
x = [real(zfrom), real(zto)];
y = [imag(zfrom), imag(zto)];

plot(x,y,plot_prop,'LineWidth',1.2*0+1)
end

