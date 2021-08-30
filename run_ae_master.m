%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           RUN AE MASTER
%
%           Created by: Erik Toller,
%                       erik.toller@geo.uu.se
%                       Department of Earth Sciences, Uppsala University
%                       SWEDEN
%
%           Last updated: 2021-08-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
close all
clc
%% ASSINING THE GLOBAL PARAMETERS
sigma_11inf = .2*0;
nu = 0.3;
G = 20e3/10000;
kappa = 3 - 4 * nu;

% load('data_files/10000_crack_data.mat')
% z1a = z1(1:end-6);
% z2a = z2(1:end-6);
% load('data_files/3_intersecting_cracks.mat')
% z1a = [z1a,z1];
% z2a = [z2a,z2];

m_IT = [500:250:1000];
for IT = 1:length(m_IT)
% if NIT == 2
%     xfrom = .1;
%     xto = .45;
%     yfrom = -.2;
%     yto = yfrom -(xfrom-xto);
%     collect = real(z1a)<xto & real(z1a)>xfrom & imag(z1a)<yto & imag(z1a)>yfrom;
%     z1a = z1a(collect);
%     z2a = z2a(collect);
% end
% if NIT == 3
%     xfrom = .25-.05;
%     xto = .35+.05;
%     yfrom = -.1-.05;
%     yto = yfrom -(xfrom-xto);
%     collect = real(z1a)<xto & real(z1a)>xfrom & imag(z1a)<yto & imag(z1a)>yfrom;
%     z1a = z1a(collect);
%     z2a = z2a(collect);
% end

%% ANALYTIC ELEMENT FOR A CRACK
z1a = [complex(.6,-.3), complex(-.5,-.5)];
z2a = [complex(-.6,.3), complex(.5,.5)];
% z1a = [complex(.5,.3), complex(0,-.3)];
% z2a = [complex(-.5,.3), complex(0,.29)];
% z1a = [complex(-.6,-.6),complex(-.6,.6)];
% z2a = [complex(.6,.6),complex(.6,-.6)];
% z1a = [complex(-.6,-.8),complex(.6,-.8)];
% z2a = [complex(.6,.8),complex(-.6,.8)];
% z1a = [complex(0,-.5),complex(.5,.25)];
% z2a = [complex(0,.5),complex(-.5,.25)];
% z1a = [complex(.2,-.5),complex(-.5,.2)];
% z2a = [complex(.2,.5),complex(.5,.2)];
% z1a = [complex(-.6,0),complex(0,-.6)];
% z2a = [complex(.6,0),complex(0,.6)];
% z1a = [complex(-.6,-.6),complex(-.6,.6*0)];
% z2a = [complex(.1,.1),complex(.1,-.1)];
% z1a = [complex(0,-.6)];
% z2a = [complex(0,.6)];
% z1a = [complex(-.5,-.5),complex(.6,-.2),complex(.6,.1)];
% z2a = [complex(.5,.5),complex(-.6,.2),complex(.1,-.4)];
% z1a = [complex(-.6,-.6),complex(-.4,.4),complex(.6,.6),complex(.4,-.4)];
% z2a = [complex(-.4,.4),complex(.6,.6),complex(.4,-.4),complex(-.6,-.6)];
% z1a = [complex(-.4,-.4),complex(-.4,.4),complex(.4,-.4),complex(.4,.4)];
% z2a = [complex(0,0),complex(0,0),complex(0,0),complex(0,0)];
% z1a = [complex(-.4,-.6),complex(-.6,.4),complex(.4,.6),complex(.6,-.4)];
% z2a = [complex(-.4,.6),complex(.6,.4),complex(.4,-.6),complex(-.6,-.4)];
% z1a = [complex(.4,-.6),complex(-.6,.4)];
% z2a = [complex(.4,.6),complex(.6,.4)];
na = length(z1a);

%% COEFFICIENTS FOR AE
% Analaytic Element
ma = m_IT(IT);
Na = ma*2*0+2000;
p = -1;

cond = 1e-15;
NIT = 300;
numa = Na;

%% PLOT DIMENSIONS AND RESOLUTION
% Dimensions
xfrom = -1;
xto = 1;
yfrom = -1;
yto = yfrom -(xfrom-xto);
% if NIT == 2
%     xfrom = .1;
%     xto = .45;
%     yfrom = -.2;
%     yto = yfrom -(xfrom-xto);
% end
% if NIT == 3
%     xfrom = .25;
%     xto = .35;
%     yfrom = -.1;
%     yto = yfrom -(xfrom-xto);
% end
% Resolution
Nx = 200;
Ny = Nx;
Nw = 50;
Ntraj = 400;
lvs_traj = 30;
lvs = 40;

%% EXPORT THE INPUT DATA
% Analytic Element
if na == 0
    z1a = [];
    z2a = [];
    pa = [];
    La = [];
    mua = [];
else
    pa = p.*ones(1,na);
    La = zeros(1,na);
    mua = zeros(1,na);
    for ii = 1:na
        La(ii) = sqrt((z2a(ii)-z1a(ii))*conj(z2a(ii)-z1a(ii)));
        mua(ii) = angle(z2a(ii)-z1a(ii));
    end
end

% Tarjectories start line
if xfrom < 0
    xtraj = [xfrom+yfrom*1.5*1i, xfrom+yto*1.5*1i];
    ytraj = [xfrom*1.5+yfrom*1i, xto*1.5+yfrom*1i];
elseif xfrom > 0
    xtraj = [xfrom+yfrom*1.5*1i, xfrom+yto*1.5*1i];
    ytraj = [xfrom/1.5+yfrom*1i, xto*1.5+yfrom*1i];
end
xtraj_vec = []; ytraj_vec = [];
for ii = 1:2
    xtraj_vec = [xtraj_vec, real(xtraj(ii)), imag(xtraj(ii))];
    ytraj_vec = [ytraj_vec, real(ytraj(ii)), imag(ytraj(ii))];
end

% Write the bin-files for C++ program
A = [sigma_11inf,kappa,G,na,ma,Na,cond, NIT,...
    real(z1a),imag(z1a),real(z2a),imag(z2a),pa,La,mua]; % Vector to write
input_file = fopen('geometry_data.bin','w');
fwrite(input_file, A, 'double');
fclose(input_file);

B = [xfrom,xto,yfrom,yto,Nx,Ny,Nw,Ntraj,lvs_traj,xtraj_vec,ytraj_vec]; % Vector to write
plot_file = fopen('plot_data.bin','w');
fwrite(plot_file, B, 'double');
fclose(plot_file);

%% RUN THE C++ PROGRAM
system('run_AE_LE_master.exe');

%% IMPORT THE RESULTS
% Load the BIN-files
data_file = fopen('data.bin', 'r');
dim_file = fopen('dim_data.bin', 'r');
coef_file = fopen('input_data.bin', 'r');
[A,~] = fread(data_file,'double');
[B,~] = fread(dim_file,'double');
[C,~] = fread(coef_file,'double');
fclose(data_file);
fclose(dim_file);
fclose(coef_file);
% Collect the results
Nx = B(1);
Ny = B(2);
Nw = B(3);
Ntraj = B(4);
lvs_traj = B(5);
na = B(6);
error_med_a_re = B(7);
error_mean_a_re = B(8);
error_max_a_re = B(9);
error_med_a_im = B(10);
error_mean_a_im = B(11);
error_max_a_im = B(12);
error_med_int = B(13);
error_mean_int = B(14);
error_max_int = B(15);
z1a = zeros(1,na);
z2a = zeros(1,na);
pos = 15;
for ii = 1:na
    re = pos + ii;
    im = pos + na + ii;
    z1a(ii) = complex(B(re),B(im));
    re = pos + na*2 + ii;
    im = pos + na*3 + ii;
    z2a(ii) = complex(B(re),B(im));
end
pos = pos + 4*na;
% Collect the grids
x_vec = A(1:Nx);
y_vec = A((Nx+1):(Nx+Ny));
x_vecw = A((Nx+Ny+1):(Nx+Ny+Nw));
y_vecw = A((Nx+Ny+Nw+1):(Nx+Ny+Nw+Nw));
grid_11 = zeros(Nx,Ny);
grid_22 = zeros(Nx,Ny);
grid_12 = zeros(Nx,Ny);
grid_1 = zeros(Nx,Ny);
grid_2 = zeros(Nx,Ny);
theta_p = zeros(Nx,Ny);
traj_1 = zeros(lvs_traj*2,Ntraj);
traj_2 = zeros(lvs_traj*2,Ntraj);
grid_w = zeros(Nw,Nw);
traj_w = zeros(Nw,Ntraj);
T_check = zeros(1,numa*na);
start = Nx+Ny+Nw+Nw+1;
for ii = 1:Nx
    stop = start + Nx -1;
    grid_11(:,ii) = A(start:stop);
    start = stop + 1;
end
for ii = 1:Nx
    stop = start + Nx -1;
    grid_22(:,ii) = A(start:stop);
    start = stop + 1;
end
for ii = 1:Nx
    stop = start + Nx -1;
    grid_12(:,ii) = A(start:stop);
    start = stop + 1;
end
for ii = 1:Nx
    stop = start + Nx -1;
    grid_1(:,ii) = A(start:stop);
    start = stop + 1;
end
for ii = 1:Nx
    stop = start + Nx -1;
    grid_2(:,ii) = A(start:stop);
    start = stop + 1;
end
for ii = 1:Nx
    stop = start + Nx -1;
    theta_p(:,ii) = A(start:stop);
    start = stop + 1;
end
for ii = 1:lvs_traj*2
    stop = start + Ntraj*2 -1;
    for jj = 1:Ntraj
        re = start + jj - 1;
        im = start + Ntraj + jj - 1;
        traj_1(ii,jj) = complex(A(re),A(im));
    end
    start = stop + 1;
end
for ii = 1:lvs_traj*2
    stop = start + Ntraj*2 -1;
    for jj = 1:Ntraj
        re = start + jj - 1;
        im = start + Ntraj + jj - 1;
        traj_2(ii,jj) = complex(A(re),A(im));
    end
    start = stop + 1;
end
for ii = 1:Nw
    stop = start + Nw*2 -1;
    for jj = 1:Nw
        re = start + jj - 1;
        im = start + Nw + jj - 1;
        grid_w(jj,ii) = complex(A(re),A(im));
    end
    start = stop + 1;
end
for ii = 1:Nw*2
    stop = start + Ntraj*2 -1;
    for jj = 1:Ntraj
        re = start + jj - 1;
        im = start + Ntraj + jj - 1;
        traj_w(ii,jj) = complex(A(re),A(im));
    end
    start = stop + 1;
end
for ii = 1:1
    stop = start + numa -1;
    for jj = 1:numa*na
        re = start + jj - 1;
        im = start + numa + jj - 1;
        T_check(jj) = complex(A(re),A(im));
    end
    start = stop + 1;
end

% Get the coefficients
sigma_11inf = C(1);
kappa = C(2);
G = C(3);
% na = C(4);
% ma = C(5);
a = zeros(na,ma);
start = 6+na*6;
for ii = 1:na
    stop = start + ma*2 -1;
    for jj = 1:ma
        re = start + jj - 1;
        im = start + ma + jj - 1;
        a(ii,jj) = complex(C(re),C(im));
    end
    start = stop + 1;
end

% Arrange the trajectories into a sinlge vector
traj_1p = [flip(traj_1(1:lvs_traj,:),2),traj_1(lvs_traj+1:end,:)];
traj_2p = [flip(traj_2(1:lvs_traj,:),2),traj_2(lvs_traj+1:end,:)];

%% PLOT THE RESULTS
% Set font to LaTeX
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',15)
% Display set contour levels
disp(' ')
disp(['Contour levels set to: ',num2str(lvs)])
disp(' ')

% Create the figures
% disp('figure ( 1/10)')
% create_figure()
% contour(x_vec, y_vec, grid_1,lvs,'blue','LineWidth',1.0);
% for ii = 1:na
%     Plot_line(z1a(ii),z2a(ii),'black')
% end
% legend('$\sigma_{1}$')
% axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])
% 
% disp('figure ( 2/10)')
% create_figure()
% contour(x_vec, y_vec, grid_2,lvs,'red','LineWidth',1.0);
% for ii = 1:na
%     Plot_line(z1a(ii),z2a(ii),'black')
% end
% legend('$\sigma_{2}$')
% axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])
% 
% disp('figure ( 3/10)')
% create_figure()
% contour(x_vec, y_vec, grid_11,lvs,'blue','LineWidth',1.0);
% for ii = 1:na
%     Plot_line(z1a(ii),z2a(ii),'black')
% end
% legend('$\sigma_{11}$')
% axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])
% 
% disp('figure ( 4/10)')
% create_figure()
% contour(x_vec, y_vec, grid_22,lvs,'red','LineWidth',1.0);
% for ii = 1:na
%     Plot_line(z1a(ii),z2a(ii),'black')
% end
% legend('$\sigma_{22}$')
% axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])
% 
% disp('figure ( 5/10)')
% create_figure()
% contour(x_vec, y_vec, grid_12,lvs,'green','LineWidth',1.0);
% for ii = 1:na
%     Plot_line(z1a(ii),z2a(ii),'black')
% end
% legend('$\sigma_{12}$')
% axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])
% 
% disp('figure ( 6/10)')
% create_figure()
% contour(x_vec, y_vec, theta_p,lvs,'red','LineWidth',1.0);
% for ii = 1:na
%     Plot_line(z1a(ii),z2a(ii),'black')
% end
% legend('$\theta_{p}$')
% axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])
% 
% disp('figure ( 7/10)')
% create_figure()
% for ii = 1:lvs_traj
%     p1 = plot(real(traj_1p(ii,:)),imag(traj_1p(ii,:)),'blue','LineWidth',1.0);
%     p2 = plot(real(traj_2p(ii,:)),imag(traj_2p(ii,:)),'red','LineWidth',1.0);
% end
% for ii = 1:na
%     Plot_line(z1a(ii),z2a(ii),'black')
% end
% legend([p1 p2], '$\sigma_{1}$','$\sigma_{2}$')
% axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])
% 
% disp('figure ( 8/10)')
% create_figure()
% quiver(x_vecw, y_vecw, real(grid_w), -imag(grid_w),'blue');
% for ii = 1:Nw*2
%     p1 = plot(real(traj_w(ii,:)),imag(traj_w(ii,:)),'blue');
% end
% for ii = 1:na
%     Plot_line(z1a(ii),z2a(ii),'black')
% end
% legend('$w$')
% axis([x_vecw(1) x_vecw(end) y_vecw(1) y_vecw(end)])
% 
% disp('figure ( 9/10)')
% create_figure()
% contour(x_vec, y_vec, grid_1,lvs,'blue','LineWidth',1.0);
% for ii = 1:lvs_traj
%     p1 = plot(real(traj_1p(ii,:)),imag(traj_1p(ii,:)),': blue','LineWidth',1.0);
% end
% for ii = 1:na
%     Plot_line(z1a(ii),z2a(ii),'black')
% end
% legend('$\sigma_{1}$')
% axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])
% 
% disp('figure (10/10)')
% create_figure()
% contour(x_vec, y_vec, grid_2,lvs,'red','LineWidth',1.0);
% for ii = 1:lvs_traj
%     p1 = plot(real(traj_2p(ii,:)),imag(traj_2p(ii,:)),': red','LineWidth',1.0);
% end
% for ii = 1:na
%     Plot_line(z1a(ii),z2a(ii),'black')
% end
% legend('$\sigma_{2}$')
% axis([x_vec(1) x_vec(end) y_vec(1) y_vec(end)])


%%
error_medim(IT) = abs(error_med_a_im)/abs(2*p);
error_medre(IT) = abs(error_med_a_re)/abs(2*p);
error_med(IT) = abs(error_med_int)/abs(2*p);
error_meandim(IT) = abs(error_mean_a_im)/abs(2*p);
error_meandre(IT) = abs(error_mean_a_re)/abs(2*p);
error_meand(IT) = abs(error_mean_int)/abs(2*p);
IT_vec(IT) = ma;

% figure
% hold on
% for ii = 1:na
%     T_temp = T_check(1+numa*(ii-1):numa*ii);
%     plot(1:round(numa),real(T_temp),'blue')
%     plot(1:round(numa),imag(T_temp),'red')
% end

%% SAVE WORKPLACE DATA
% Save the data as todays date + version
today = date;
str_file = ['data_files\simulation_',today,'_v0.mat'];
while isfile(str_file)
    vpos = find(str_file == 'v');
    version = str2double(str_file((vpos(end)+1):end-4))+1;
    str_file = [str_file(1:vpos(end)),num2str(version),'.mat'];
end
save(str_file,...
    'G','grid_1','grid_11','grid_12','grid_2','grid_22','grid_w',...
    'kappa','La','lvs','lvs_traj','ma','mua',...
    'na','Na','Ntraj','nu','Nw','Nx','Ny','pa',...
    'sigma_11inf','theta_p','traj_1','traj_1p','traj_2','traj_2p',...
    'traj_w','x_vec','x_vecw','xfrom','xto','xtraj','xtraj_vec',...
    'y_vec','y_vecw','yfrom','yto','ytraj','ytraj_vec','z1a','z2a',...
    'error_max_a_im','error_max_a_re',...
    'error_mean_a_im','error_mean_a_re',...
    'error_med_a_im','error_med_a_re',...
    'error_mean_int','error_med_int','error_max_int')

disp(' ')
disp(['Data saved as: ','"',str_file,'"'])
end
figure
semilogy(IT_vec, (error_medre),'blue',IT_vec, (error_medim),'red',...
    IT_vec, (error_med),'black',IT_vec, (error_meandre),'blue-.',...
    IT_vec, (error_meandim),'red-.',IT_vec, (error_meand),'black-.')
%legend('$|\tau^{11}|$','$|t_n-p|$','$|t_s|$')
grid minor

figure
loglog(IT_vec, (error_medre),'blue',IT_vec, (error_medim),'red',...
    IT_vec, (error_med),'black',IT_vec, (error_meandre),'blue-.',...
    IT_vec, (error_meandim),'red-.',IT_vec, (error_meand),'black-.')
%legend('$|\tau^{11}|$','$|t_n-p|$','$|t_s|$')
grid minor



