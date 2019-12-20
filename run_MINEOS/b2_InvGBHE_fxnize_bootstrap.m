% Code to invert for (Gc, Gs); (Bc, Bs); (Hc, Hs); (Ec, Es) in a 1D profile
% as in Russell et al.; JGR (2019). The methodology and kernel scaling 
% directly follows Montagner & Nataf 1986.
%
%plot native
clear;
close all;

% Define full dataset (fds) values as global
global data_fds kernels_fds 
%%
%%%%%%% INVERSION PARAMETERS %%%%%%%
% Define model depth
par.model_depth_G = 400;
par.model_depth_B = par.model_depth_G;
par.model_depth_H = par.model_depth_B;
par.model_depth_E = 25; %25;

% Ignore crust...
par.Vpv_cutoff = 7800; % 4000 (crust) [km/s] cut kernels where velocities lower than this value

% smoothing (second derivative)
    par.alphF_GBH = 4e2; %9e2; 
    par.alphF_E = 6e2; %3e3; 
    
% Flatness (first derivative);
    par.alphJ_GBH = 1e1;
    par.alphJ_E = 2e2; %2e3; 
    
% Damp below certain depth
    par.G_DAMP = 5e0; %5e1;
    par.dep_zero_damp = 300; %300 %[km]

% Norm Damping
    par.alphH_GBH = 5e0; %1e0;  %smaller fit more
    par.alphH_E = 5e0;  %smaller fit more
    
% Scaling Ratios
    par.epsilonBG = 1e5; %1e3; %1e0;  % Enforce B/G ratio
    par.BGratio = 1.25; %1.25; %1.25;
    par.epsilonHG = 1e5; %1e3; %1e0;  % Enforce H/G ratio
    par.HGratio = -0.11; %0.11; %0.25; %0.25; %1.25;
    
% General data weighting
par.eps_d = 1e0; % data weighting
par.eps_f = 1e0; % constraint weighting

% Break into layers
par.is_brk = 0; % Break model into layers?
par.z_brks = [25 100 300]; % [km] depths to form breaks
par.alphF_brk = 2e2; %first derivative smoothing
par.alphJ_brk = 2e2; %second derivative smoothing
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrap parameters
par.nbs = 1; %1000 % Number of bootstrap iterations
par.nbins = 60; % for heat plots
par.chi2red_thresh = 1.5; % chi2 threshold; will only consider models with chi2 less than this value when determining confidence intervals

isfig = 0; % save figures?
is_RMS = 0; % Use RMS uncertainties for anisotropy strength measurements rather than proper 95 confidence?
issavemat = 0; % Save results .mat file?

ylims = [0 300];
par.APM = 114; % absolute plate motion (GSRM 2.1; NNR) https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion-calculator.html
par.FSD = 75; % fossil spreading direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load General Parameters
% parameter_ACFLN
parameter_FRECHET
DATAmat_path = param.DATAmat_path;
modelpath = param.modelpath;
if par.nbs==1
    par.chi2red_thresh = 9999;
end

%% Initialize model structure
Velmodel = load(modelpath);
I = find(Velmodel.card.vpv > par.Vpv_cutoff);
iz_cut = length(Velmodel.card.vpv)-I(end);

model = model_init(Velmodel.card,par,iz_cut);

% Make directories
if par.Vpv_cutoff == 1500
    top = 'noh20';
elseif par.Vpv_cutoff == 3000
    top = 'noh20sed';
elseif par.Vpv_cutoff>6000 && par.Vpv_cutoff<9000
    top = 'noh20sedcrust';
else
    top = 'test';
end
figpath = [param.figpath,top,'_dep',num2str(par.model_depth_G),'/'];
if exist(figpath,'dir')==0
    mkdir(figpath);
end

%% Do bootstrap

% [ data_bs ] = index_dataset(data, par.nbs);
[ data_bs ] = bootstrap_dataset(data, par.nbs);

% Begin bootstrap iterations
model_save = model;
for ibs = 1:par.nbs
    model = model_save;
    
    %% Get Mineos dispersion
    model = read_Mineos_disp(data_bs(ibs),model,param);
%     if ibs == 1
%         model_fds = model;
%     end

    %% Convert (A, phi) anisotropy measurements to (Acos, Asin)
    data_bs(ibs) = Aphi2sincos(data_bs(ibs),is_RMS);
    if ibs == 1
        data_fds = data_bs(ibs);
    end
    %% Setup kernels for inversion
    kernels = load_kernels(param, par, data_bs(ibs), model);
    if ibs == 1
        kernels_fds = kernels;
    end
    
    %% Do the inversions

    % Rayleigh & Love 2-theta inversions
    model = run_GBHinv( data_bs(ibs), model, kernels, par, ibs );

    % Love 4-theta inversions
    model = run_Einv( data_bs(ibs), model, kernels, par, ibs );
    
    model_bs(ibs) = model;
    kernels_bs(ibs) = kernels;
    
    if mod(ibs,50)==0
        disp(ibs);
    end
end
for ibs = 1:par.nbs
    model_bs(ibs).G.G_dir(model_bs(ibs).G.G_dir<30) = model_bs(ibs).G.G_dir(model_bs(ibs).G.G_dir<30) + 180;
end


[ensemble] = collect_models(model_bs, data_bs(1), par, {'G','B','H','E'});

[stats, ensemble] = get_stats(ensemble, par);

%% Get Indices for plotting

I_R0 = data.rayl.mode_br_ani == 0;
I_R1 = data.rayl.mode_br_ani == 1;
I_L0 = data.love.mode_br_ani == 0;
I_L1 = data.love.mode_br_ani == 1;

A2_L = [data.love.A2(I_L0), data.love.A2(I_L1)];
A4_L = [data.love.A4(I_L0), data.love.A4(I_L1)];
err_A2_L = [data.love.err_2A(I_L0), data.love.err_2A(I_L1)];
err_A4_L = [data.love.err_4A(I_L0), data.love.err_4A(I_L1)];
phi2_L = [data.love.phi2(I_L0), data.love.phi2(I_L1)];
phi4_L = [data.love.phi4(I_L0), data.love.phi4(I_L1)];
err_phi2_L = [data.love.err_phi2(I_L0), data.love.err_phi2(I_L1)];
err_phi4_L = [data.love.err_phi4(I_L0), data.love.err_phi4(I_L1)];
Tperiods = [data.love.periods_ani(I_L0), data.love.periods_ani(I_L1)];

%% Plot Histograms
if par.nbs > 1
f53 = figure(53);
set(gcf,'position',[27           1        1185         664]);

[G_cnt, G_bins] = plot_hist([2,2,1],ensemble.X2red_2,par.nbs);
title('G','fontsize',18);
% [G_cnt_bsiter, G_bins_bsiter] = plot_hist([2,2,2],chi2_GBH_red_bsiter,nbs);
% title('G (bootstraps)','fontsize',18);
[E_cnt, E_bins] = plot_hist([2,2,3],ensemble.X2red_4,par.nbs);
title('E','fontsize',18);
% [E_cnt_bsiter, E_bins_bsiter] = plot_hist([2,2,4],chi2_C_red_bsiter,nbs);
% title('E (bootstraps)','fontsize',18);

%% Heat Maps 1
figure(51); clf;
set(gcf,'position',[1           1        1250         697]);
max_cbar = 0.95;
ncont = 10;
LBLFNT=18
cmapp = parula;

%G
subplot(1,4,1); box on;
contourf(stats.G.G_mag_bins*100,model.G.z,stats.G.G_mag_cnts'./max(stats.G.G_mag_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(stats.G.G_mag_l95*100,model.G.z,'-y','linewidth',3);
plot(stats.G.G_mag_u95*100,model.G.z,'-y','linewidth',3);
plot(stats.G.G_mag_l68*100,model.G.z,'-r','linewidth',3);
plot(stats.G.G_mag_u68*100,model.G.z,'-r','linewidth',3);
ylim([12 300]);
caxis([0 max_cbar]);
ylabel('Depth (km)','fontsize',18)
xlabel('Strength (%)','fontsize',18)
title('G/L','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');
colormap(cmapp);

subplot(1,4,2); box on;
contourf(stats.G.G_dir_bins,model.G.z,stats.G.G_dir_cnts'./max(stats.G.G_dir_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(stats.G.G_dir_l95,model.G.z,'-y','linewidth',3);
plot(stats.G.G_dir_u95,model.G.z,'-y','linewidth',3);
plot(stats.G.G_dir_l68,model.G.z,'-r','linewidth',3);
plot(stats.G.G_dir_u68,model.G.z,'-r','linewidth',3);
% ylabel('Depth (km)','fontsize',18)
xlabel('Azimuth (\circ)','fontsize',18)
title('G/L','fontsize',18);
ylim([12 300]);
xlim([50 200]);
caxis([0 max_cbar]);
plot([ 78 78],[0 250],'--','color',[0.8 0.8 0.8],'Linewidth',3);
APM = 296.74;
plot([ APM-180 APM-180 ],[0 400],'--','color',[1 1 1], 'Linewidth',3);
set(gca,'YDir','reverse')
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');
set(gca,'xtick',[0:30:180]);

% C
subplot(1,4,3); box on;
contourf(stats.E.E_mag_bins*100,model.E.z,stats.E.E_mag_cnts'./max(stats.E.E_mag_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(stats.E.E_mag_l95*100,model.E.z,'-y','linewidth',2);
plot(stats.E.E_mag_u95*100,model.E.z,'-y','linewidth',2);
plot(stats.E.E_mag_l68*100,model.E.z,'-r','linewidth',2);
plot(stats.E.E_mag_u68*100,model.E.z,'-r','linewidth',2);
% ylim([12 40]);
caxis([0 max_cbar]);
ylim([12 300])
xlim([min(stats.G.G_mag_bins*100), max(stats.G.G_mag_bins*100)]);
xlabel('Strength (%)','fontsize',18)
title('E/N','fontsize',18);
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');

subplot(1,4,4); box on;
contourf(stats.E.E_dir_bins,model.E.z,stats.E.E_dir_cnts'./max(stats.E.E_dir_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(stats.E.E_dir_l95,model.E.z,'-y','linewidth',2);
plot(stats.E.E_dir_u95,model.E.z,'-y','linewidth',2);
plot(stats.E.E_dir_l68,model.E.z,'-r','linewidth',2);
plot(stats.E.E_dir_u68,model.E.z,'-r','linewidth',2);
% ylabel('Depth (km)','fontsize',18)
xlabel('Azimuth (\circ)','fontsize',18)
title('E/N','fontsize',18);
caxis([0 max_cbar]);
% ylim([12 40]);
ylim([12 300]);
xlim([50 150]);
plot([ 78 78],[0 250],'--','color',[0.8 0.8 0.8],'Linewidth',3);
APM = 296.74;
plot([ APM-180 APM-180 ],[0 400],'--','color',[1 1 1], 'Linewidth',3);
set(gca,'YDir','reverse')
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');
set(gca,'xtick',[0:30:180]);

% %%
end

%% plot inversion result for P2P strength and fastdir as function of Depth (SPLIT)
figure(103);clf;
reds = brewermap(5,'reds');
blues = brewermap(5,'blues');
clr_G = reds(4,:);
clr_G_bs = reds(2,:);
clr_E = blues(4,:);
clr_E_bs = blues(2,:);

% STRENGTH
Nx = 2; Ny = 1;
sidegap = 0.10; topgap = 0.10; botgap = 0.10; vgap = 0.05; hgap = 0.05; cbar_bot = 0.04;
width = (1 - vgap*(Nx-1)-2*sidegap)/Nx; height = (1 - topgap - botgap - (Ny-1)*hgap)/Ny;
set(gcf,'position',[ 118   290   1000   750]);
% set(gcf,'position',[210     1   830   704]);
set(gcf,'color','w');
LBLFNT=18
ix = 1; iy = 1;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height]);
hold on; set(gca,'fontsize',18,'linewidth',2);
% h1(2) = plot(model.B.B_mag*100,model.B.z,'-','color',[1 0.7 0.7],'linewidth',4);hold on;
% h1(3) = plot(model.H.H_mag*100,model.H.z,'-','color',[0.7 0.7 1],'linewidth',4);hold on;
plot(ensemble.G.G_mag(:,ensemble.I_good_X2red2)*100,model.G.z,'-','color',clr_G_bs,'linewidth',2);hold on;
plot(ensemble.E.E_mag(:,ensemble.I_good_X2red4)*100,model.E.z,'-','color',clr_E_bs,'linewidth',2);hold on;
h1(1) = plot(stats.G.G_mag_med*100,model.G.z,'-','color',clr_G,'linewidth',4);hold on;
h1(4) = plot(stats.E.E_mag_med*100,model.E.z,'-','color',clr_E,'linewidth',4);hold on;
xlim([0 6]); 
ylim([0 300]);
ylabel('Depth (km)','fontsize',18)
xlabel('Strength (%)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'YDir','reverse')
% legend(h1,{'G/L (2\theta)','B/A (2\theta)','H/F (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off','Location','southeast');



% AZIMUTH
ix = 1;iy = 2;
left = sidegap + (iy-1)*(vgap+width);bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height]);
hold on;set(gca,'fontsize',18,'linewidth',2);hold on;set(gca,'fontsize',18,'linewidth',2);
% h1(2) = plot(model.B.B_dir, model.B.z,'-','color',[1 0.7 0.7],'linewidth',4);hold on;
% h1(3) = plot(model.H.H_dir, model.H.z,'-','color',[0.7 0.7 1],'linewidth',4);hold on;
plot(ensemble.G.G_dir(:,ensemble.I_good_X2red2),model.G.z,'-','color',clr_G_bs,'linewidth',2);hold on;
plot(ensemble.E.E_dir(:,ensemble.I_good_X2red4),model.E.z,'-','color',clr_E_bs,'linewidth',2);hold on;
h1(1) = plot(stats.G.G_dir_med,model.G.z,'-','color',clr_G,'linewidth',4);hold on;
h1(4) = plot(stats.E.E_dir_med,model.E.z,'-','color',clr_E,'linewidth',4);hold on;
% plot(stats.G.G_dir_bs(:,1),model.G.z,'-k','linewidth',2);
plot([ 78 78],[0 250],'k--','Linewidth',3);
plot([ par.APM par.APM ],[50 400],'--','color',[0.6 0.6 0.6], 'Linewidth',3);
set(gca,'YDir','reverse')
% xlim([20 200]);
ylim([0 300]);
xlabel('Azimuth (\circ)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
set(gca,'xtick',[0:30:180]);

if isfig
    save2pdf([figpath,'azi_depth_GBHE.pdf'],103,1000)
%     export_fig([figpath,'azi_depth_GBHE_2.pdf'],'-pdf','-q100','-p0.02','-painters',103)
end

%% plot A and Phi fits (horizontal)
figure(105);clf;

% STRENGTH
set(gcf,'position',[210     1   830   704],'color','w');
subplot(2,1,1);
hold on; set(gca,'fontsize',18,'linewidth',2);
plot(model_bs(1).R0.periods,ensemble.R0_fds.A2(:,ensemble.I_good_X2red2)*200,'-','color',clr_G_bs,'linewidth',2);
plot(model_bs(1).R1.periods-.1,ensemble.R1_fds.A2(:,ensemble.I_good_X2red2)*200,'-','color',clr_G_bs,'linewidth',2);
plot(model_bs(1).L.periods,ensemble.L_fds.A2(:,ensemble.I_good_X2red2)*200,'-','color',clr_G_bs,'linewidth',2);
plot(model_bs(1).L.periods,ensemble.L_fds.A4(:,ensemble.I_good_X2red4)*200,'-','color',clr_E_bs,'linewidth',2);
plot(data.rayl.periods_ani(I_R0),data.rayl.A2(I_R0)*2*100,'o','color',clr_G,'linewidth',2,'markersize',12,'markerfacecolor',clr_G); hold on;
plot(data.rayl.periods_ani(I_R1)-.1,data.rayl.A2(I_R1)*2*100,'o','color',clr_G,'linewidth',2,'markersize',12,'markerfacecolor',clr_G); hold on;
plot(Tperiods,A2_L*2*100,'^','color',clr_G,'linewidth',2,'markersize',12,'markerfacecolor',clr_G);
plot(Tperiods,A4_L*2*100,'^','color',clr_E,'linewidth',2,'markersize',12,'markerfacecolor',clr_E);
h = errorbar(data.rayl.periods_ani(I_R0),data.rayl.A2(I_R0)*2*100,data.rayl.err_2A(I_R0)*100,'vertical','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(data.rayl.periods_ani(I_R1)-.1,data.rayl.A2(I_R1)*2*100,data.rayl.err_2A(I_R1)*100,'vertical','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(Tperiods,A2_L*2*100,err_A2_L*100,'vertical','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(Tperiods,A4_L*2*100,err_A4_L*100,'vertical','.','linewidth',2,'color','k');
h.CapSize = 0;
plot(model_bs(1).R0.periods,stats.R0.A2_med*2*100,'-','color',clr_G,'linewidth',4);
plot(model_bs(1).R1.periods-.1,stats.R1.A2_med*2*100,'-','color',clr_G,'linewidth',4);
plot(model_bs(1).L.periods,stats.L.A2_med*2*100,'-','color',clr_G,'linewidth',4);
plot(model_bs(1).L.periods,stats.L.A4_med*2*100,'-','color',clr_E,'linewidth',4);
% xlabel('Period (s)','fontsize',12);
ylabel('2A (%)','fontsize',12);
set(gca,'fontsize',12,'linewidth',1.5)
xlim([min(data.rayl.periods_ani(I_R1))-1 max(data.rayl.periods_ani(I_R0))+10])
ylim([0 5]);
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'XScale','log')


% AZIMUTH
subplot(2,1,2);
hold on;set(gca,'fontsize',18,'linewidth',2);hold on;set(gca,'fontsize',18,'linewidth',2);
plot(model_bs(1).R0.periods,ensemble.R0_fds.phi2(:,ensemble.I_good_X2red2),'-','color',clr_G_bs,'linewidth',2);
plot(model_bs(1).R1.periods-.1,ensemble.R1_fds.phi2(:,ensemble.I_good_X2red2),'-','color',clr_G_bs,'linewidth',2);
plot(model_bs(1).L.periods,ensemble.L_fds.phi2(:,ensemble.I_good_X2red2),'-','color',clr_G_bs,'linewidth',2);
plot(model_bs(1).L.periods,ensemble.L_fds.phi4(:,ensemble.I_good_X2red4),'-','color',clr_E_bs,'linewidth',2);
plot(data.rayl.periods_ani(I_R0),data.rayl.phi2(I_R0),'o','color',clr_G,'linewidth',2,'markersize',12,'markerfacecolor',clr_G);hold on;
plot(data.rayl.periods_ani(I_R1),data.rayl.phi2(I_R1),'o','color',clr_G,'linewidth',2,'markersize',12,'markerfacecolor',clr_G);hold on;
plot(Tperiods,phi2_L,'^','color',clr_G,'linewidth',2,'markersize',12,'markerfacecolor',clr_G);hold on;
plot(Tperiods,phi4_L,'^','color',clr_E,'linewidth',2,'markersize',12,'markerfacecolor',clr_E);hold on;
h = errorbar(data.rayl.periods_ani(I_R0),data.rayl.phi2(I_R0),data.rayl.err_phi2(I_R0),'vertical','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(data.rayl.periods_ani(I_R1),data.rayl.phi2(I_R1),data.rayl.err_phi2(I_R1),'vertical','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(Tperiods,phi2_L,err_phi2_L,'vertical','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(Tperiods,phi4_L,err_phi4_L,'vertical','.','linewidth',2,'color','k');
h.CapSize = 0;
plot(model_bs(1).R0.periods,stats.R0.phi2_med,'-','color',clr_G,'linewidth',4);
plot(model_bs(1).R1.periods-.1,stats.R1.phi2_med,'-','color',clr_G,'linewidth',4);
plot(model_bs(1).L.periods,stats.L.phi2_med,'-','color',clr_G,'linewidth',4);
plot(model_bs(1).L.periods,stats.L.phi4_med,'-','color',clr_E,'linewidth',4);
plot([20 400],[ par.APM par.APM ],'--','color',[0.6 0.6 0.6], 'Linewidth',3);
plot([1 250],[ par.FSD par.FSD],'k--','Linewidth',3);
plot([1 10],[ par.FSD par.FSD]+45,'k--','Linewidth',3);
plot([1 10],[ par.FSD par.FSD]+90,'k--','Linewidth',3);
xlim([min(data.rayl.periods_ani(I_R1))-1 max(data.rayl.periods_ani(I_R0))+10])
ylim([20 200]);
xlabel('Period (s)','fontsize',12);
ylabel('\psi (\circ)','fontsize',12); hold on;
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
set(gca,'xtick',[0:30:210],'XScale','log');
% legend(h1,{'G/L (2\theta)','B/A (2\theta)','H/F (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');
% legend(h1,{'G/L (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');

if isfig
    save2pdf([figpath,'BOOT_azi_depth_GBHE_bootstrap_alphF',num2str(par.alphF_GBH,'%2.1e'),'_nbs',num2str(nbs),'_ORALS_lonpar.G_DAMP0_NORM_perc_MN86_A_PHI_errprop_patch_2.pdf'],104,1000)
%     export_fig([figpath,'BOOT_azi_depth_GBHE_bootstrap_alphF',num2str(par.alphF_GBH,'%2.1e'),'_nbs',num2str(nbs),'_ORALS_lonpar.G_DAMP0_NORM_perc_MN86_A_PHI_errprop_patch2.pdf'],'-pdf','-q100','-p0.02','-painters',104)
end

%% Save results
if issavemat
    matpath = ['./mats_bootstrap/'];
    if ~exist(matpath)
        mkdir(matpath);
    end
    param.par = par;
    results.data = data;
    results.kernels = kernels_bs(1);
    results.kernels_bs = kernels_bs;
    results.model = model_bs(1);
    results.model_bs = model_bs;
    results.param = param;
    results.ensemble = ensemble;
    results.stats = stats;
    save([matpath,param.CARDID,'_nbs',num2str(par.nbs),'.mat'],'results');
end
