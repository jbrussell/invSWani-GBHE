% Plot ACFLN kernels for S0 and S1 from 5-145 s
%
% JBR - 4/28/17
%

clear all;
close all;

%% Parameters
SONLY_vec = [1]; %1 spheroidal (0 toroidal)
TONLY_vec = [0]; %1 Toroidal (0 spheroidal)
branch_vec = [0]; %0 fundamental, 1 first overtone
fid = fopen('mode.txt','w');
fprintf(fid,'%.0f %.0f %.0f',SONLY_vec,TONLY_vec,branch_vec);
fclose(fid);

parameter_FRECHET

FRECHETPATH = param.frechetpath;
DATAmat_path = param.DATAmat_path;
Gmat_S0path = param.Gmat_S0path;
Gmat_S1path = param.Gmat_S1path;
Gmat_T0path = param.Gmat_T0path;
S0periods = flip(param.S0periods);
S1periods = flip(param.S1periods);
T0periods = flip(param.T0periods);
modelpath = param.modelpath;
% periods = param.periods;
figpath = param.figpath;

ylims = [0 100]; %[11 300];
ylims_shall = [11 24];

model_depth = 400;


%% get G matrix
T0_Gmat = load(Gmat_T0path);
S0_Gmat = load(Gmat_S0path);
S1_Gmat = load(Gmat_S1path);

%% Load Velocity Model
Velmodel = load(modelpath);

%qeqeqe = 5; % no invert water;
%qeqeqe = 7; %no invert water + sediments; %7%5% no inverst for water
% qeqeqe = 12; %10; %12; %13 % no invert for water + sediments + crust;
vpv_cutoff = 8000;
I = find(Velmodel.card.vpv > vpv_cutoff);
qeqeqe = length(Velmodel.card.vpv)-I(end);

rad = Velmodel.card.rad;
Idep = find(rad >= rad(end)-model_depth*1000);
Idep = Idep(1:end-qeqeqe);
Idep2 = Idep;
rad = flipud(rad(Idep));
Vsv = flipud(Velmodel.card.vsv(Idep));%/1000);
Vsh = flipud(Velmodel.card.vsh(Idep));%/1000);
Vpv = flipud(Velmodel.card.vpv(Idep));%/1000);
Vph = flipud(Velmodel.card.vph(Idep));%/1000);
eta = flipud(Velmodel.card.eta(Idep));%/1000);
rho = flipud(Velmodel.card.rho(Idep));%/1000);

%% Get MINEOS dispersion
CARDID = param.CARDID;
qname_S = [CARDID,'.s0to200.q'];
qname_T = [CARDID,'.t0to200.q'];

%S0
mode = 0;
% Load phase and group velocities
qfile = [FRECHETPATH,qname_S];
[phV,grV,phVq] = readMINEOS_qfile2(qfile,S0periods,mode);
c_S0 = flip(phV);
U_S0 = flip(grV);

%S1
mode = 1;
% Load phase and group velocities
qfile = [FRECHETPATH,qname_S];
[phV,grV,phVq] = readMINEOS_qfile2(qfile,S1periods,mode);
c_S1 = flip(phV);
U_S1 = flip(grV);

%T0
mode = 0;
% Load phase and group velocities
qfile = [FRECHETPATH,qname_T];
[phV,grV,phVq] = readMINEOS_qfile2(qfile,T0periods,mode);
c_T0 = flip(phV);
U_T0 = flip(grV);

%% set up kernels

Idep = find(S0_Gmat.Gmatrix.rad >= S0_Gmat.Gmatrix.rad(1)-model_depth*1000);
Idep = Idep(qeqeqe+1:end);

z = (S0_Gmat.Gmatrix.rad(1)-S0_Gmat.Gmatrix.rad(Idep))/1000;

% 2 theta
% G (L) Rayleigh, (-L) Love
GGG2_R0 = S0_Gmat.Gmatrix.L(:,Idep).*(c_S0'./U_S0'.*rho'.*Vsv'.^2);
GGG2_R1 = S1_Gmat.Gmatrix.L(:,Idep).*(c_S1'./U_S1'.*rho'.*Vsv'.^2);
GGG2_L = -T0_Gmat.Gmatrix.L(:,Idep).*(c_T0'./U_T0'.*rho'.*Vsv'.^2);

% B (A) Rayleigh
BBB2_R0 = S0_Gmat.Gmatrix.A(:,Idep).*(c_S0'./U_S0'.*rho'.*Vph'.^2);
BBB2_R1 = S1_Gmat.Gmatrix.A(:,Idep).*(c_S1'./U_S1'.*rho'.*Vph'.^2);

% 4 theta
% C (A) Rayleigh, (-N) Love
CCC4 = [%S_Gmat.Gmatrix.A(:,IdepC);
    -T0_Gmat.Gmatrix.N(:,Idep)].*(c_T0'./U_T0'.*rho'.*Vsh'.^2);


%% Plot Kernels
clr = [1 0 0; 0 1 0; 0 0 1; 0.5 0.5 0.5];
clr = lines(5);

% Long Periods
fig1 = figure(1); clf;
% set(gcf,'position',[210     1   485   704]);
% set(gcf,'position',[659   212   393   493]);
set(gcf,'position',[221   230   733   471]);
% set(gcf,'color','w');
LBLFNT=18
hold on; set(gca,'fontsize',18,'linewidth',2);

subplot(1,2,2);
% S0
% L
h1 = plot(S0_Gmat.Gmatrix.L(:,Idep).*(c_S0(:)./U_S0(:).*rho'.*Vsv'.^2),z,'linewidth',3,'color',clr(1,:)); hold on;
% A
h2 = plot(S0_Gmat.Gmatrix.A(:,Idep).*(c_S0(:)./U_S0(:).*rho'.*Vph'.^2),z,'linewidth',3,'color',clr(2,:)); hold on;
% F
h3 = plot(S0_Gmat.Gmatrix.F(:,Idep).*(c_S0(:)./U_S0(:).*rho'.*eta'.*(Vph'.^2 - 2*Vsv'.^2)),z,'linewidth',3,'color',clr(3,:)); hold on;

h = [h1(1), h2(1), h3(1)];
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
legend(h,{'S0 K''_L (2\theta)','S0 K''_A (2\theta)','S0 K''_F (2\theta)'},'location','southeast','fontsize',15,'box','off');

ylabel('Depth (km)','fontsize',LBLFNT);
title('15-150 s');
ylim([11 300]);
%         xlim([1e-16 1e-11]);
set(gca,'fontsize',LBLFNT,'linewidth',1.5,'ydir','reverse')
% xlim([0 5])
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'YDir','reverse')


subplot(1,2,1)
% Short Periods
% set(gcf,'position',[210     1   485   704]);
% set(gcf,'position',[659   212   393   493]);
% set(gcf,'color','w');
LBLFNT=18
hold on; set(gca,'fontsize',18,'linewidth',2);

% S1
% L
h1 = plot(S1_Gmat.Gmatrix.L(:,Idep).*(c_S1(:)./U_S1(:).*rho'.*Vsv'.^2),z,'linewidth',3,'color',clr(1,:)); hold on;
% A
h2 = plot(S1_Gmat.Gmatrix.A(:,Idep).*(c_S1(:)./U_S1(:).*rho'.*Vph'.^2),z,'linewidth',3,'color',clr(2,:)); hold on;
% F
h3 = plot(S1_Gmat.Gmatrix.F(:,Idep).*(c_S1(:)./U_S1(:).*rho'.*eta'.*(Vph'.^2 - 2*Vsv'.^2)),z,'linewidth',3,'color',clr(3,:)); hold on;

% T0
% L
h4 = plot(T0_Gmat.Gmatrix.L(:,Idep).*(c_T0(:)./U_T0(:).*rho'.*Vsv'.^2),z,'--','linewidth',3,'color',clr(4,:)); hold on;
% N
h5 = plot(T0_Gmat.Gmatrix.N(:,Idep).*(c_T0(:)./U_T0(:).*rho'.*Vsh'.^2),z,'--','linewidth',3,'color',clr(5,:)); hold on;

h = [h1(1), h2(1), h3(1), h4(1), h5(1)];
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
legend(h,{'S1 K''_L (2\theta)','S1 K''_A (2\theta)','S1 K''_F (2\theta)','T0 K''_L (2\theta)','T0 K''_N (4\theta)'},'location','southeast','fontsize',15,'box','off');
% ylabel('Depth (km)','fontsize',15);
ylim([11 40]);
title('5-7.5 s');
%         xlim([1e-16 1e-11]);
ylabel('Depth (km)','fontsize',LBLFNT);
set(gca,'fontsize',LBLFNT,'linewidth',1.5,'ydir','reverse')
% xlim([0 5])
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'YDir','reverse')

CARDID = param.CARDID;
TYPEID = param.TTYPEID;
%print('-painters','-dpdf','-r400',[EIGPATH,CARDID,'.',TYPEID,'.',num2str(j),'mod.',num2str(N_modes),'_fix.pdf']);
save2pdf([FRECHETPATH,'ACFLNkernels_scaled_S0S1T0_GBHE_ORALS','.pdf'],fig1,1000);
