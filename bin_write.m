function [] = bin_write(A, B)
% BIN_WRITE Wrties the binary files for C++ computations
%   This program provides the program stress_field_calculator.cpp with the
%   necessary bin-files. It can either be called with the data collected in
%   two vectors or it can be run by itself.
%
%   The length of the vectors are coded in and needs to be updated when
%   more elements are added.

% Assigning variables
sigma_11inf = 5;
nu = 0.3;
G = 20e3;
kappa = 3 - 4 * nu;

load('200_crack_data.mat')
m = 20;
nc = length(z1);
% p = -1.*rand(1,nc);
L = zeros(1,nc);
mu = zeros(1,nc);
for ii = 1:nc
    L(ii) = sqrt((z2(ii)-z1(ii))*conj(z2(ii)-z1(ii)));
    mu(ii) = angle(z2(ii)-z1(ii));
end

N = m*2;

% Plot dim and resolution
% xfrom = 0;
% xto = 0.25;
% yfrom = -0.05;
% yto = 0.20;
xfrom = -1;
xto = 1;
yfrom = -1;
yto = 1;
xx = .5*(xto+xfrom);
yy = .5*(yto+yfrom);
Nx = 4;
Ny = 4;
Nw = 40;
Ntraj = 600;
lvs_traj = 80;
xtraj = [xx+yfrom*1i,xx+yto*1i];
ytraj = [xfrom+yy*1i, xto+yy*1i];
xtraj_vec = []; ytraj_vec = [];
for ii = 1:2
    xtraj_vec = [xtraj_vec, real(xtraj(ii)), imag(xtraj(ii))];
    ytraj_vec = [ytraj_vec, real(ytraj(ii)), imag(ytraj(ii))];
end

% Write the bin-files for C++ program
A = [sigma_11inf,kappa,G,nc,m,N,real(z1),imag(z1),real(z2),imag(z2),p,L,mu]; % Vector to write
input_file = fopen('geometry_data.bin','w');
fwrite(input_file, A, 'double');
fclose(input_file);

B = [xfrom,xto,yfrom,yto,Nx,Ny,Nw,Ntraj,lvs_traj,xtraj_vec,ytraj_vec]; % Vector to write
plot_file = fopen('plot_data.bin','w');
fwrite(plot_file, B, 'double');
fclose(plot_file);

disp('The output files has been written.')



end
