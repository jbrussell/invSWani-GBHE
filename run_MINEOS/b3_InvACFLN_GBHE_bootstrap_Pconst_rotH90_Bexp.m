clear;
close all;
%% inverse for Gc and Gs..
% by pylin.patty 2014.11
%
% 2/19/17 JBR
% Edit to invert for Gc Gs Cc and Cs for both Rayleigh and Love
%
% 1/12/18 JBR : ORALS
%   - In this version, I re-multiply the errors by the amount I shrank them
%   to put them back to the original scale
%
% 4/22/18 JBR : JGR
%   - New version of bootstrap resampling where balanced resampling is
%   used.
%   - This version bootstraps 23 G data (18 S0-2theta + 5 S1-2theta + 5 T0-2theta
%                              5 E data (5 T0-4theta) 
%
% 4/30/19 JBR : 
%   - Realized that H should be 90º from G,B,E. This is easy to impose by
%   simply making the H/G constraint < 0. It doesn't actually affect the
%   answer much since the ratio is so small, but makes more sense because
%   the H kernel is negative and now it doesn't run away with B.
%% Parameters
parameter_FRECHET

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
FRECHETPATH = param.frechetpath;

ylims = [0 300];
ylims_shall = [11 40];

%% LOAD Cij
cij = load_cij;

% Calculate parameters
cij_calc = cij_calculations(cij);

cij_ind = [2 5]; %[1 2 3 5];

%% get G matrix

T0_Gmat = load(Gmat_T0path);
S0_Gmat = load(Gmat_S0path);
S1_Gmat = load(Gmat_S1path);

% Choose card and periods to replace kernels in case they are affected by nonlinearity
if 0
    start_card = 'Nomelt_taper_eta_crust'; % Name of card to replace with
    Iper_replace = 5; % 7.5 s (index of period to replace)
    T0_Gmat_start = load([path2runMINEOS,'/MODE/FRECHET/',start_card,'/','Gmatrix_T0_',num2str(param.T0periods(1)),'_',num2str(param.T0periods(end)),'s_',start_card,'.mat']);
    T0_Gmat.Gmatrix.L(5,:) = T0_Gmat_start.Gmatrix.L(5,:);
    T0_Gmat.Gmatrix.N(5,:) = T0_Gmat_start.Gmatrix.N(5,:);
end

%model_depth = T_Gmat.Gmatrix.model_depth;
model_depth_G = 400;
model_depth_B = 400;
model_depth_H = model_depth_B;
model_depth_C = 40; %25;

dep_zero_damp = 300;
%%%%%%% INVERSION PARAMETERS %%%%%%%

% smoothing (second derivative)
    alphF_GBH = 12e2; %9e2; %2e3 
    alphF_C = 3e3; %9e2; %3e3; 
    
% Flatness (first derivative);
    alphJ_C = 2e2; %2e3; %3e3; 
    
% Damp below certain depth
    G_DAMP = 5e1; %5e1; %3e2; 

% Norm Damping
    alphH_GBH = 1e0; %0.005 %good : 0.001 %smaller fit more
    alphH_C = 0; %5e1; %0.005 %good : 0.001 %smaller fit more
    
% Scaling Ratios
    epsilonBG = 7e2; %7e2; %1e3; %1e0;  % Enforce B/G ratio
    BGratio = 1.5; %1.25; % empirical scaling ~1.2911
    epsilonHG = 1e3; %1e3; %1e0;  % Enforce H/G ratio
    HGratio = -0.11; %0.25; %0.25; %1.25; %%%% SHOULD BE NEGATIVE TO MAKE 90º ROTATED FORM G,B,E (4/30/19)
    is_downweightshallow_BG = 1;% B/G ratio; downweight shallow layers?
    
% Bootstrap
    nbs = 2000; %200; %100; % number of iterations
    nbins = 60; % for heat plots
    chi2red_thresh = 2; %1.25; % Calculate confidence bounds using only models with reduced chi^2 less than this value (< 1 overfititng; >1 underfitting)
%     strength_bins = [0:(.07-0)/100:.07]: % for heat plots
%     fastdir_bins = [50:(150-50)/100:50]: % for heat plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P Constraints at Moho
DAMP_B_A_moho_7km = 1e3; %1e3;
DAMP_E_N_moho_7km = 1e2; %3e2;
% c11 = 207.15;
% c22 = 232.71;
% c16 = 0.53;
% c26 = 2.47;
% c12_2c66 = 213.39;
% 
% B_A_7km = 7.7/100;


isfig = 0;
is_RMS = 1;

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


%% get model info for pho and Vs
Velmodel = load(modelpath);

%qeqeqe = 5; % no invert water;
%qeqeqe = 7; %no invert water + sediments; %7%5% no inverst for water
% qeqeqe = 12; %10; %12; %13 % no invert for water + sediments + crust;
vpv_cutoff = 7500; % 4000 (crust)
I = find(Velmodel.card.vpv > vpv_cutoff);
qeqeqe = length(Velmodel.card.vpv)-I(end);
if vpv_cutoff == 1500
    top = 'noh20';
elseif vpv_cutoff == 3000
    top = 'noh20sed';
elseif vpv_cutoff == 7500
    top = 'noh20sedcrust';
else
    top = 'test';
end
figpath2 = [figpath,top,'_dep',num2str(model_depth_G),'/'];
figpath = figpath2;
if exist(figpath,'dir')==0
    mkdir(figpath);
end

rad = Velmodel.card.rad;
Idep = find(rad >= rad(end)-model_depth_G*1000);
Idep = Idep(1:end-qeqeqe);
IdepB = find(rad >= rad(end)-model_depth_B*1000);
IdepB = IdepB(1:end-qeqeqe);
IdepH = IdepB;
IdepC = find(rad >= rad(end)-model_depth_C*1000);
IdepC = IdepC(1:end-qeqeqe);
rad = flipud(rad(Idep));
depthlayer_G = midpts(flipud(Velmodel.card.z(Idep))')';
depthlayer_B = midpts(flipud(Velmodel.card.z(IdepB))')';
depthlayer_H = midpts(flipud(Velmodel.card.z(IdepH))')';
depthlayer_C = midpts(flipud(Velmodel.card.z(IdepC))')';
Vsv = midpts(flipud(Velmodel.card.vsv(Idep))')';%/1000);
Vsh = midpts(flipud(Velmodel.card.vsh(Idep))')';%/1000);
Vsh_C = midpts(flipud(Velmodel.card.vsh(IdepC))')';%/1000);
Vpv = midpts(flipud(Velmodel.card.vpv(Idep))')';%/1000);
Vph = midpts(flipud(Velmodel.card.vph(Idep))')';%/1000);
eta = midpts(flipud(Velmodel.card.eta(Idep))')';%/1000);
rho = midpts(flipud(Velmodel.card.rho(Idep))')';%/1000);
rho_C = midpts(flipud(Velmodel.card.rho(IdepC))')';%/1000);

L_norm = rho .* Vsv.^2;
N_norm = rho_C .* Vsh_C.^2; 
C_norm = rho .* Vpv.^2;
A_norm = rho .* Vph.^2;
F_norm = rho .* eta.* (Vph.^2 - 2*Vsv.^2);

%% Load P-wave constraints
fid = fopen('Pconstraints.txt');
header = fgetl(fid);
scan = textscan(fid,'%s %f %f');
fclose(fid);
Pcons_moho = scan{2};
Pcons_7km = scan{3};

%% Calculate P-wave constraints
c11_moho = Pcons_moho(1);          c11_7km = Pcons_7km(1);
c22_moho = Pcons_moho(2);          c22_7km = Pcons_7km(2);
c16_moho = Pcons_moho(3);          c16_7km = Pcons_7km(3);
c26_moho = Pcons_moho(4);          c26_7km = Pcons_7km(4);
c12_2c66_moho = Pcons_moho(5);     c12_2c66_7km = Pcons_7km(5);
Vph_moho = mean(Pcons_moho(6:7));  Vph_7km = mean(Pcons_7km(6:7));
%%%%%%% Constraints on B/A %%%%%%%
% Moho B
A_moho = 3/8*(c11_moho+c22_moho) + 1/4*c12_2c66_moho;
rho_moho = A_moho/Vph_moho^2 * 1000;
Bc_A_moho = 1/2*(c11_moho-c22_moho) /A_moho;
Bs_A_moho = (c16_moho+c26_moho) /A_moho;
B_A_moho = sqrt(Bc_A_moho^2+Bs_A_moho^2);
phi_B_moho = 0.5*atan2d(Bs_A_moho,Bc_A_moho);

% 7km B
A_7km = 3/8*(c11_7km+c22_7km) + 1/4*c12_2c66_7km;
rho_7km = A_7km/Vph_7km^2 * 1000;
Bc_A_7km = 1/2*(c11_7km-c22_7km) /A_7km;
Bs_A_7km = (c16_7km+c26_7km) /A_7km;
B_A_7km = sqrt(Bc_A_7km^2+Bs_A_7km^2);
phi_B_7km = 0.5*atan2d(Bs_A_7km,Bc_A_7km);

% Interpolate from Moho to moho+7km
I_moho_7km = find(depthlayer_B<depthlayer_B(1)+7);
depth_moho_7km = depthlayer_B(I_moho_7km);
Bc_A_moho_7km = linspace(Bc_A_moho,Bc_A_7km,length(depth_moho_7km));
Bs_A_moho_7km = linspace(Bs_A_moho,Bs_A_7km,length(depth_moho_7km));
B_A_moho_7km = linspace(B_A_moho,B_A_7km,length(depth_moho_7km));
phi_B_moho_7km = linspace(phi_B_moho,phi_B_7km,length(depth_moho_7km));
% rho_7km = rho(I_moho_7km(end));

%%%%%%% Constraints on E/N %%%%%%%
% Moho E
N_moho = N_norm(1);
Ec_N_moho = (1/8*(c11_moho+c22_moho) - 1/4*c12_2c66_moho) /N_moho * 1e9; % 0.0218
Es_N_moho = 1/2*(c16_moho-c26_moho) /N_moho * 1e9; % -0.0130
E_N_moho = sqrt(Ec_N_moho^2+Es_N_moho^2);
phi_E_moho = 0.25*atan2d(Es_N_moho,Ec_N_moho);

% 7km E
N_7km = N_norm(I_moho_7km(end));
Ec_N_7km = (1/8*(c11_7km+c22_7km) - 1/4*c12_2c66_7km) /N_7km * 1e9; % 0.0218
Es_N_7km = 1/2*(c16_7km-c26_7km) /N_7km * 1e9; % -0.0130
E_N_7km = sqrt(Ec_N_7km^2+Es_N_7km^2);
phi_E_7km = 0.25*atan2d(Es_N_7km,Ec_N_7km);

% Interpolate from Moho to moho+7km
Ec_N_moho_7km = linspace(Ec_N_moho,Ec_N_7km,length(depth_moho_7km));
Es_N_moho_7km = linspace(Es_N_moho,Es_N_7km,length(depth_moho_7km));
E_N_moho_7km = linspace(E_N_moho,E_N_7km,length(depth_moho_7km));
phi_E_moho_7km = linspace(phi_E_moho,phi_E_7km,length(depth_moho_7km));

%% get the p2p and theta for each period from aniso.mat

load(DATAmat_path);
aniso_R = data.rayl;
aniso_L = data.love;
all_periods = aniso_R.periods;

% Rayleigh
c_iso_R = aniso_R.c_iso*1000;
phi2_R = aniso_R.phi2;
phi4_R = aniso_R.phi4;
err_phi2_R = aniso_R.err_phi2;
err_phi4_R = aniso_R.err_phi4;
A2_R_perc = aniso_R.A2;
A4_R_perc = aniso_R.A4;
err_A2_R_perc = aniso_R.err_2A;
err_A4_R_perc = aniso_R.err_4A;
A2_R = aniso_R.A2;%.* c_iso_R;
A4_R = aniso_R.A4;%.* c_iso_R;
err_A2_R = aniso_R.err_2A;%.* c_iso_R;
err_A4_R = aniso_R.err_4A;%.* c_iso_R;

% Love
c_iso_L = aniso_L.c_iso*1000;
phi2_L = aniso_L.phi2;
phi4_L = aniso_L.phi4;
err_phi2_L = aniso_L.err_phi2;
err_phi4_L = aniso_L.err_phi4;
A2_L_perc = aniso_L.A2;
A4_L_perc = aniso_L.A4;
err_A2_L_perc = aniso_L.err_2A;
err_A4_L_perc = aniso_L.err_4A;
A2_L = aniso_L.A2;%.* c_iso_L;
A4_L = aniso_L.A4;%.* c_iso_L;
err_A2_L = aniso_L.err_2A;%.* c_iso_L;
err_A4_L = aniso_L.err_4A;%.* c_iso_L;

err_pct = 0.6;
if is_RMS
    err_A2_R_perc = aniso_R.wRMS_2A/2 * err_pct;
    err_A4_R_perc = aniso_R.wRMS_4A/2 * err_pct;
    err_A2_L_perc = aniso_L.wRMS_2A/2 * err_pct;
    err_A4_L_perc = aniso_L.wRMS_4A/2 * err_pct;
    
    err_A2_R = aniso_R.wRMS_2A/2 * err_pct;%.*c_iso_R;
    err_A4_R = aniso_R.wRMS_4A/2 * err_pct;%.*c_iso_R;
    err_A2_L = aniso_L.wRMS_2A/2 * err_pct;%.*c_iso_L;
    err_A4_L = aniso_L.wRMS_4A/2 * err_pct;%.*c_iso_L;
end
err_A2_R(12:13) = err_A2_R(12:13) * 3.5;
err_A2_R_perc(12:13) = err_A2_R_perc(12:13) * 3.5;
err_phi2_R(12:13) = err_phi2_R(12:13) * 3.5;
% err_phi2_R(18) = err_phi2_R(18) * 3.5; % 5 seconds



% Convert to Montagner and Nataf form
% Rayleigh
c2_R = A2_R.*cosd(2*phi2_R); c2_R_save = c2_R;
c4_R = A4_R.*cosd(4*phi4_R); c4_R_save = c4_R;
s2_R = A2_R.*sind(2*phi2_R); s2_R_save = s2_R;
s4_R = A4_R.*sind(4*phi4_R); s4_R_save = s4_R;
c2_std_R = sqrt( (err_A2_R.*cosd(2*phi2_R)).^2 + (-err_phi2_R*pi/180*2.*A2_R.*sind(2*phi2_R)).^2 ); c2_std_R_save = c2_std_R;
c4_std_R = sqrt( (err_A4_R.*cosd(4*phi4_R)).^2 + (-err_phi4_R*pi/180*4.*A4_R.*sind(4*phi4_R)).^2 ); c4_std_R_save = c4_std_R;
s2_std_R = sqrt( (err_A2_R.*sind(2*phi2_R)).^2 + ( err_phi2_R*pi/180*2.*A2_R.*cosd(2*phi2_R)).^2 ); s2_std_R_save = s2_std_R;
s4_std_R = sqrt( (err_A4_R.*sind(4*phi4_R)).^2 + ( err_phi4_R*pi/180*4.*A4_R.*cosd(4*phi4_R)).^2 ); s4_std_R_save = s4_std_R;

% Love
c2_L = A2_L.*cosd(2*phi2_L); c2_L_save = c2_L;
c4_L = A4_L.*cosd(4*phi4_L); c4_L_save = c4_L;
s2_L = A2_L.*sind(2*phi2_L); s2_L_save = s2_L;
s4_L = A4_L.*sind(4*phi4_L); s4_L_save = s4_L;
c2_std_L = sqrt( (err_A2_L.*cosd(2*phi2_L)).^2 + (-err_phi2_L*pi/180*2.*A2_L.*sind(2*phi2_L)).^2 ); c2_std_L_save = c2_std_L;
c4_std_L = sqrt( (err_A4_L.*cosd(4*phi4_L)).^2 + (-err_phi4_L*pi/180*4.*A4_L.*sind(4*phi4_L)).^2 ); c4_std_L_save = c4_std_L;
s2_std_L = sqrt( (err_A2_L.*sind(2*phi2_L)).^2 + ( err_phi2_L*pi/180*2.*A2_L.*cosd(2*phi2_L)).^2 ); s2_std_L_save = s2_std_L;
s4_std_L = sqrt( (err_A4_L.*sind(4*phi4_L)).^2 + ( err_phi4_L*pi/180*4.*A4_L.*cosd(4*phi4_L)).^2 ); s4_std_L_save = s4_std_L;

%% Do bootstrap
ndata_S0 = length(S0periods);
ndata_S1 = length(S1periods);
ndata_S = ndata_S0 + ndata_S1;
ndata_T = length(T0periods);

% Setup indexes for balanced random resampling
[ indxs_ST ] = balanced_resampling( ndata_S + ndata_T,nbs-1 );
indxs_ST = sort(indxs_ST,1);
indxs_ST = [[1:ndata_S + ndata_T]' , indxs_ST];
[ indxs_T ] = balanced_resampling( ndata_T,nbs-1 );
indxs_T = sort(indxs_T,1);
indxs_T = [[1:ndata_T]' , indxs_T];
for ibs = 1:nbs
% I_bs_S = bootstrapindex(ndata_S,nperm_S);
I_bs_S = indxs_ST(indxs_ST(:,ibs)<=ndata_S,ibs); % index for S0 and S1
I_bs_flip_S = sort(ndata_S+1-I_bs_S); % index for flipped periods
I_bs_S0 = I_bs_S(I_bs_S <= length(S0periods));
I_bs_flip_S0 = sort(ndata_S0+1-I_bs_S0); % index for flipped periods
I_bs_S1 = I_bs_S(I_bs_S > length(S0periods)) - ndata_S0;
I_bs_flip_S1 = sort(ndata_S1+1-I_bs_S1); % index for flipped periods
% I_bs_T = bootstrapindex(ndata_T,nperm_T);
% I_bs_T = I_bs_S1;
I_bs_T2thet = indxs_ST(indxs_ST(:,ibs)>ndata_S,ibs)-ndata_S; % index for T-2theta 
I_bs_flip_T2thet = sort(ndata_T+1-I_bs_T2thet);
I_bs_T4thet = indxs_T(:,ibs);
I_bs_flip_T4thet = sort(ndata_T+1-I_bs_T4thet); % index for flipped periods

%% Make G matrix for inversion

I_R0 = aniso_R.mode_br == 0;
I_R1 = aniso_R.mode_br == 1;
c_iso_R0 = c_iso_R(I_R0);

IdepG = find(S0_Gmat.Gmatrix.rad >= S0_Gmat.Gmatrix.rad(1)-model_depth_G*1000);
IdepG = IdepG(qeqeqe+1:end);
IdepB = find(S0_Gmat.Gmatrix.rad >= S0_Gmat.Gmatrix.rad(1)-model_depth_B*1000);
IdepB = IdepB(qeqeqe+1:end);
IdepC = find(T0_Gmat.Gmatrix.rad >= T0_Gmat.Gmatrix.rad(1)-model_depth_C*1000);
IdepC = IdepC(qeqeqe+1:end);
IdepH = IdepB;

radG = S0_Gmat.Gmatrix.rad(IdepG);
% drG = -flip([flip(diff(radG))]);
drG = -[diff(radG)];
drB = drG;
drH = drB;

radC = T0_Gmat.Gmatrix.rad(IdepC);
drC = -[diff(radC)];

% 2 theta
% G (L) Rayleigh, (-L) Love
GGG2_R0 = midpts(S0_Gmat.Gmatrix.L(I_bs_flip_S0,IdepG)).*drG'.*(c_S0(I_bs_flip_S0)'./U_S0(I_bs_flip_S0)'.*rho'.*Vsv'.^2);
GGG2_R1 = midpts(S1_Gmat.Gmatrix.L(I_bs_flip_S1,IdepG)).*drG'.*(c_S1(I_bs_flip_S1)'./U_S1(I_bs_flip_S1)'.*rho'.*Vsv'.^2);
GGG2_L = midpts(-T0_Gmat.Gmatrix.L(I_bs_flip_T2thet,IdepG)).*drG'.*(c_T0(I_bs_flip_T2thet)'./U_T0(I_bs_flip_T2thet)'.*rho'.*Vsv'.^2);
GGG2_R0_allper = midpts(S0_Gmat.Gmatrix.L(:,IdepG)).*drG'.*(c_S0'./U_S0'.*rho'.*Vsv'.^2);
GGG2_R1_allper = midpts(S1_Gmat.Gmatrix.L(:,IdepG)).*drG'.*(c_S1'./U_S1'.*rho'.*Vsv'.^2);
GGG2_L_allper = midpts(-T0_Gmat.Gmatrix.L(:,IdepG)).*drG'.*(c_T0'./U_T0'.*rho'.*Vsv'.^2);

% B (A) Rayleigh
BBB2_R0 = midpts(S0_Gmat.Gmatrix.A(I_bs_flip_S0,IdepB)).*drB'.*(c_S0(I_bs_flip_S0)'./U_S0(I_bs_flip_S0)'.*rho'.*Vph'.^2);
BBB2_R1 = midpts(S1_Gmat.Gmatrix.A(I_bs_flip_S1,IdepB)).*drB'.*(c_S1(I_bs_flip_S1)'./U_S1(I_bs_flip_S1)'.*rho'.*Vph'.^2);
BBB2_R0_allper = midpts(S0_Gmat.Gmatrix.A(:,IdepB)).*drB'.*(c_S0'./U_S0'.*rho'.*Vph'.^2);
BBB2_R1_allper = midpts(S1_Gmat.Gmatrix.A(:,IdepB)).*drB'.*(c_S1'./U_S1'.*rho'.*Vph'.^2);

% H (F) Rayleigh
HHH2_R0 = midpts(S0_Gmat.Gmatrix.F(I_bs_flip_S0,IdepH)).*drH'.*(c_S0(I_bs_flip_S0)'./U_S0(I_bs_flip_S0)'.*rho'.*eta'.*(Vph'.^2 - 2*Vsv'.^2));
HHH2_R1 = midpts(S1_Gmat.Gmatrix.F(I_bs_flip_S1,IdepH)).*drH'.*(c_S1(I_bs_flip_S1)'./U_S1(I_bs_flip_S1)'.*rho'.*eta'.*(Vph'.^2 - 2*Vsv'.^2));
HHH2_R0_allper = midpts(S0_Gmat.Gmatrix.F(:,IdepH)).*drH'.*(c_S0'./U_S0'.*rho'.*eta'.*(Vph'.^2 - 2*Vsv'.^2));
HHH2_R1_allper = midpts(S1_Gmat.Gmatrix.F(:,IdepH)).*drH'.*(c_S1'./U_S1'.*rho'.*eta'.*(Vph'.^2 - 2*Vsv'.^2));


% 4 theta
% C (A) Rayleigh, (-N) Love
CCC4 = [%S_Gmat.Gmatrix.A(:,IdepC);
    midpts(-T0_Gmat.Gmatrix.N(I_bs_flip_T4thet,IdepC))].*drC'.*(c_T0(I_bs_flip_T4thet)'./U_T0(I_bs_flip_T4thet)'.*rho_C'.*Vsh_C'.^2);
CCC4_allper = [%S_Gmat.Gmatrix.A(:,IdepC);
    midpts(-T0_Gmat.Gmatrix.N(:,IdepC))].*drC'.*(c_T0'./U_T0'.*rho_C'.*Vsh_C'.^2);

[~, nlayer_G ] = size(midpts(S0_Gmat.Gmatrix.L(:,IdepG)));
nperiod_G = size(S0_Gmat.Gmatrix.L(:,IdepG),1)+size(S1_Gmat.Gmatrix.L(:,IdepG),1);
nperiods_G = nperiod_G;
nlayers_G = nlayer_G;
[~, nlayer_B ] = size(midpts(S0_Gmat.Gmatrix.A(:,IdepB)));
nperiod_B = size(S0_Gmat.Gmatrix.A(:,IdepB),1)+size(S1_Gmat.Gmatrix.A(:,IdepB),1);
nperiods_B = nperiod_B;
nlayers_B = nlayer_B;
[nperiod_C, nlayer_C ] = size(midpts(T0_Gmat.Gmatrix.N(:,IdepC)));
nperiods_C = nperiod_C;
nlayers_C = nlayer_C;
nlayers_H = nlayers_B;

if S0periods(1)>S0periods(end)
    GGG2_R0 = flip(GGG2_R0,1);
    GGG2_R1 = flip(GGG2_R1,1);
    GGG2_L = flip(GGG2_L,1);
    BBB2_R0 = flip(BBB2_R0,1);
    BBB2_R1 = flip(BBB2_R1,1);
    HHH2_R0 = flip(HHH2_R0,1);
    HHH2_R1 = flip(HHH2_R1,1);
    CCC4 = flip(CCC4,1);
    
    GGG2_R0_allper = flip(GGG2_R0_allper,1);
    GGG2_R1_allper = flip(GGG2_R1_allper,1);
    GGG2_L_allper = flip(GGG2_L_allper,1);
    BBB2_R0_allper = flip(BBB2_R0_allper,1);
    BBB2_R1_allper = flip(BBB2_R1_allper,1);
    HHH2_R0_allper = flip(HHH2_R0_allper,1);
    HHH2_R1_allper = flip(HHH2_R1_allper,1);
    CCC4_allper = flip(CCC4_allper,1);
end
GGG2 = [GGG2_R0;
        GGG2_R1;
        GGG2_L];

%% Build data vectors
I_R0 = aniso_R.mode_br == 0;
I_R1 = aniso_R.mode_br == 1;

% 2 theta - Cosine part
c2_R0 = c2_R_save(I_R0);
c2_R1 = c2_R_save(I_R1);
    c2_R = [c2_R0(I_bs_S0), c2_R1(I_bs_S1)]';
c2_L = c2_L_save(I_bs_T2thet)';
c2_std_R0 = c2_std_R_save(I_R0);
c2_std_R1 = c2_std_R_save(I_R1);
    c2_std_R = [c2_std_R0(I_bs_S0), c2_std_R1(I_bs_S1)]';
c2_std_L = c2_std_L_save(I_bs_T2thet)';
% 2 theta - Sine part
s2_R0 = s2_R_save(I_R0);
s2_R1 = s2_R_save(I_R1);
    s2_R = [s2_R0(I_bs_S0), s2_R1(I_bs_S1)]';
s2_L = s2_L_save(I_bs_T2thet)';
s2_std_R0 = s2_std_R_save(I_R0);
s2_std_R1 = s2_std_R_save(I_R1);
    s2_std_R = [s2_std_R0(I_bs_S0), s2_std_R1(I_bs_S1)]';
s2_std_L = s2_std_L_save(I_bs_T2thet)';
% 4 theta - Cosine part
c4_R = c4_R_save';
c4_L = c4_L_save(I_bs_T4thet)';
c4_std_R = c4_std_R_save';
c4_std_L = c4_std_L_save(I_bs_T4thet)';
% 4 theta - Sine part
s4_R = s4_R_save';
s4_L = s4_L_save(I_bs_T4thet)';
s4_std_R = s4_std_R_save';
s4_std_L = s4_std_L_save(I_bs_T4thet)';

%% INVERSION

% Setup damping to zero beneath specified depth (JOSH)
GBH_EYE = eye(nlayers_G+nlayers_B+nlayers_H);
I_DAMP = find(midpts(S0_Gmat.depthlayer(IdepG)')' >= dep_zero_damp);
G_EYE_DAMP = [GBH_EYE(I_DAMP,:)];

nlayer_G_DAMP = size(G_EYE_DAMP,1);

% Damp linearly from 0 to G_DAMP
G_DAMP_vec = linspace(0,G_DAMP,nlayer_G_DAMP)';

% P-wave constraints on B/A
nlayersPconstr = length(Bc_A_moho_7km);
B_A_Pconstr = eye(nlayers_B);
B_A_Pconstr = B_A_Pconstr(I_moho_7km,:);
%% G,B,H (RAYLEIGH)

% damping matrix
H00 = eye(3*nlayer_G);
% smoothing matrix
F00 = 2*eye(1*nlayer_G);
Fup = -1*[zeros(1*nlayer_G-1,1) eye(1*nlayer_G-1) ;zeros(1,1*nlayer_G) ];
Fdown = -1*[zeros(1,1*nlayer_G); eye(1*nlayer_G-1) zeros(1*nlayer_G-1,1) ];
F00 = F00+Fup+Fdown;
F00(1,:) = 0; F00(end,:) = 0;
F00_G = [F00, zeros(nlayer_G), zeros(nlayer_G)];
F00_B = [zeros(nlayer_G), F00, zeros(nlayer_G)];
F00_H = [zeros(nlayer_G), zeros(nlayer_G), F00];

if is_downweightshallow_BG
    std_BG = ones(nlayers_G,1)./(linspace(0.5,1,nlayers_G).^(3)'*epsilonBG); % B/G ratio; downweight shallow layers
else
    std_BG = ones(nlayers_G,1)/epsilonBG; % B/G ratio
end

% Cosine
d_GBH_c2 = [c2_R;
            c2_L;
           zeros(nlayer_G_DAMP,1);
           zeros(nlayers_B,1); % B/G ratio
           zeros(nlayers_H,1); % H/G ratio
           Bc_A_moho_7km'; % B/A P-wave constraint from moho to moho+7km
           ];

GGGBBBHHH2 = [GGG2_R0, BBB2_R0, HHH2_R0;
              GGG2_R1, BBB2_R1, HHH2_R1;
              GGG2_L, zeros(size(GGG2_L)), zeros(size(GGG2_L));
            G_EYE_DAMP;
           -BGratio*eye(nlayers_G), eye(nlayers_B), zeros(nlayers_H); % B/G ratio
           -HGratio*eye(nlayers_G), zeros(nlayers_B), eye(nlayers_H); % H/G ratio
           zeros(size(B_A_Pconstr))     , B_A_Pconstr     , zeros(size(B_A_Pconstr)); % B/A P-wave constraint from moho to moho+7km
           ];

% invert for MGc
std_GBH_c2 = [c2_std_R;
              c2_std_L;
             ones(nlayer_G_DAMP,1)./G_DAMP_vec;
%              ones(nlayers_G,1)/epsilonBG; % B/G ratio
%              ones(nlayers_G,1)./(linspace(0.5,1,nlayers_G).^(3)'*epsilonBG); % B/G ratio; downweight shallow layers
             std_BG;
             ones(nlayers_G,1)/epsilonHG; % H/G ratio
             ones(nlayersPconstr,1)/DAMP_B_A_moho_7km; % B/A P-wave constraint from moho to moho+7km
             ];
C_std_GBH_c2 = diag(std_GBH_c2.^2) ; % std matrix
H0_GBH = H00./norm(H00).*norm(GGGBBBHHH2');
F0_G = F00_G./norm(F00_G).*norm(GGGBBBHHH2');
F0_B = F00_B./norm(F00_B).*norm(GGGBBBHHH2');
F0_H = F00_H./norm(F00_H).*norm(GGGBBBHHH2');
H_GBH = alphH_GBH * H0_GBH;
F_G = alphF_GBH * F0_G;
F_B = alphF_GBH * F0_B;
F_H = alphF_GBH * F0_H;
% Q = H'*H + F'*F;
Q_GBH = F_G'*F_G + F_B'*F_B + F_H'*F_H + H_GBH'*H_GBH;
GGBBHH2 = inv(GGGBBBHHH2'* inv(C_std_GBH_c2)*GGGBBBHHH2+ Q_GBH) * GGGBBBHHH2';
MGBH_c2 = GGBBHH2*inv(C_std_GBH_c2)*d_GBH_c2;

MG_c2 = MGBH_c2(1:nlayers_G);
MB_c2 = MGBH_c2(nlayers_G+1:nlayers_G+nlayers_B);
MH_c2 = MGBH_c2(nlayers_G+nlayers_B+1:nlayers_G+nlayers_B+nlayers_H);
% display(MB_c2./MG_c2);
% display(MH_c2./MG_c2);

new_c2_G_R0 = GGG2_R0_allper*MG_c2;
new_c2_B_R0 = BBB2_R0_allper*MB_c2;
new_c2_H_R0 = HHH2_R0_allper*MH_c2;
new_c2_GBH_R0 = new_c2_G_R0 + new_c2_B_R0 + new_c2_H_R0;
new_c2_GBH_R0_test = GGGBBBHHH2*MGBH_c2;
res_c2_GBH_R0 = c2_R0' - new_c2_GBH_R0;
new_c2_G_R0_bsiter = GGG2_R0*MG_c2;
new_c2_B_R0_bsiter = BBB2_R0*MB_c2;
new_c2_H_R0_bsiter = HHH2_R0*MH_c2;
new_c2_GBH_R0_bsiter = new_c2_G_R0_bsiter + new_c2_B_R0_bsiter + new_c2_H_R0_bsiter;
res_c2_GBH_R0_bsiter = c2_R0(I_bs_S0)' - new_c2_GBH_R0_bsiter;

new_c2_G_R1 = GGG2_R1_allper*MG_c2;
new_c2_B_R1 = BBB2_R1_allper*MB_c2;
new_c2_H_R1 = HHH2_R1_allper*MH_c2;
new_c2_GBH_R1 = new_c2_G_R1 + new_c2_B_R1 + new_c2_H_R1;
new_c2_GBH_R1_test = GGGBBBHHH2*MGBH_c2;
res_c2_GBH_R1 = c2_R1' - new_c2_GBH_R1;
new_c2_G_R1_bsiter = GGG2_R1*MG_c2;
new_c2_B_R1_bsiter = BBB2_R1*MB_c2;
new_c2_H_R1_bsiter = HHH2_R1*MH_c2;
new_c2_GBH_R1_bsiter = new_c2_G_R1_bsiter + new_c2_B_R1_bsiter + new_c2_H_R1_bsiter;
res_c2_GBH_R1_bsiter = c2_R1(I_bs_S1)' - new_c2_GBH_R1_bsiter;


new_c2_G_L = GGG2_L_allper*MG_c2;
res_c2_G_L = c2_L_save' - new_c2_G_L;
new_c2_G_L_bsiter = GGG2_L*MG_c2;
res_c2_G_L_bsiter = c2_L - new_c2_G_L_bsiter;




% Sine
d_GBH_s2 = [s2_R;
            s2_L;
           zeros(nlayer_G_DAMP,1);
           zeros(nlayers_B,1); % B/G ratio
           zeros(nlayers_H,1); % H/G ratio
           Bs_A_moho_7km'; % B/A P-wave constraint from moho to moho+7km
           ];

% invert for MGs
std_GBH_s2 = [s2_std_R;
              s2_std_L;
             ones(nlayer_G_DAMP,1)./G_DAMP_vec;
%              ones(nlayers_G,1)/epsilonBG; % B/G ratio
%              ones(nlayers_G,1)./(linspace(0.5,1,nlayers_G).^(3)'*epsilonBG); % B/G ratio; downweight shallow layers
             std_BG;
             ones(nlayers_G,1)/epsilonHG; % H/G ratio
             ones(nlayersPconstr,1)/DAMP_B_A_moho_7km; % B/A P-wave constraint from moho to moho+7km
             ];
C_std_GBH_s2 = diag(std_GBH_s2.^2) ; % std matrix
GGBBHH2 = inv(GGGBBBHHH2'* inv(C_std_GBH_s2)*GGGBBBHHH2+ Q_GBH) * GGGBBBHHH2';
MGBH_s2 = GGBBHH2*inv(C_std_GBH_s2)*d_GBH_s2;

MG_s2 = MGBH_s2(1:nlayers_G);
MB_s2 = MGBH_s2(nlayers_G+1:nlayers_G+nlayers_B);
MH_s2 = MGBH_s2(nlayers_G+nlayers_B+1:nlayers_G+nlayers_B+nlayers_H);
% display(MB_s2./MG_s2);
% display(MH_s2./MG_s2);

new_s2_G_R0 = GGG2_R0_allper*MG_s2;
new_s2_B_R0 = BBB2_R0_allper*MB_s2;
new_s2_H_R0 = HHH2_R0_allper*MH_s2;
new_s2_GBH_R0 = new_s2_G_R0 + new_s2_B_R0 + new_s2_H_R0;
new_s2_GBH_R0_test = GGGBBBHHH2*MGBH_s2;
res_s2_GBH_R0 = s2_R0' - new_s2_GBH_R0;
new_s2_G_R0_bsiter = GGG2_R0*MG_s2;
new_s2_B_R0_bsiter = BBB2_R0*MB_s2;
new_s2_H_R0_bsiter = HHH2_R0*MH_s2;
new_s2_GBH_R0_bsiter = new_s2_G_R0_bsiter + new_s2_B_R0_bsiter + new_s2_H_R0_bsiter;
res_s2_GBH_R0_bsiter = s2_R0(I_bs_S0)' - new_s2_GBH_R0_bsiter;

new_s2_G_R1 = GGG2_R1_allper*MG_s2;
new_s2_B_R1 = BBB2_R1_allper*MB_s2;
new_s2_H_R1 = HHH2_R1_allper*MH_s2;
new_s2_GBH_R1 = new_s2_G_R1 + new_s2_B_R1 + new_s2_H_R1;
new_s2_GBH_R1_test = GGGBBBHHH2*MGBH_s2;
res_s2_GBH_R1 = s2_R1' - new_s2_GBH_R1;
new_s2_G_R1_bsiter = GGG2_R1*MG_s2;
new_s2_B_R1_bsiter = BBB2_R1*MB_s2;
new_s2_H_R1_bsiter = HHH2_R1*MH_s2;
new_s2_GBH_R1_bsiter = new_s2_G_R1_bsiter + new_s2_B_R1_bsiter + new_s2_H_R1_bsiter;
res_s2_GBH_R1_bsiter = s2_R1(I_bs_S1)' - new_s2_GBH_R1_bsiter;

new_s2_G_L = GGG2_L_allper*MG_s2;
res_s2_G_L = s2_L_save' - new_s2_G_L;
new_s2_G_L_bsiter = GGG2_L*MG_s2;
res_s2_G_L_bsiter = s2_L - new_s2_G_L_bsiter;

%% C
% P-wave constraints on E/N
nlayersPconstr_E = length(Ec_N_moho_7km);
E_N_Pconstr = eye(nlayers_C);
E_N_Pconstr = E_N_Pconstr(I_moho_7km,:);

% damping matrix
H0 = eye(nlayer_C);
% smoothing matrix
F0 = 2*eye(nlayer_C);
Fup = -1*[zeros(nlayer_C-1,1) eye(nlayer_C-1) ;zeros(1,nlayer_C) ];
Fdown = -1*[zeros(1,nlayer_C); eye(nlayer_C-1) zeros(nlayer_C-1,1) ];
F0 = F0+Fup+Fdown;
F0(1,:) = 0; F0(end,:) = 0;
F0(1,1) = 1; F0(1,2) = -1;
F0(end,end-1) = 1; F0(end,end) = -1;
% Flatness matrix
J0 = 1*eye(nlayer_C);
Jdown = -1*[zeros(1,nlayer_C); eye(nlayer_C-1) zeros(nlayer_C-1,1) ];
J0 = J0+Jdown;
J0(1,:) = 0; J0(end,:) = 0;
% J0(1,1) = 1; F0(1,2) = -1;
% J0(end,end-1) = 1; F0(end,end) = -1;

% Cosine
d_C_c4 = [c4_L;
          Ec_N_moho_7km';
          ];
      
CCC4_P = [CCC4;
          E_N_Pconstr;
          ];

% invert for MCc
std_C_c4 = [c4_std_L;
            ones(nlayersPconstr_E,1)/DAMP_E_N_moho_7km;
            ];
C_std_C_c4 = diag(std_C_c4.^2) ; % std matrix
H0_C = H0./norm(H0).*norm(CCC4_P');
F0_C = F0./norm(F0).*norm(CCC4_P');
J0_C = J0./norm(J0).*norm(CCC4_P');
H_C = alphH_C * H0_C;
F_C = alphF_C * F0_C;
J_C = alphJ_C * J0_C;
Q_C = F_C'*F_C + H_C'*H_C + J_C'*J_C;
CC4 = inv(CCC4_P'* inv(C_std_C_c4)*CCC4_P+ Q_C) * CCC4_P';
MC_c4 = CC4*inv(C_std_C_c4)*d_C_c4;

new_c4_C_L = CCC4_allper*MC_c4;
res_c4_C_L = c4_L_save' - new_c4_C_L;
new_c4_C_L_bsiter = CCC4*MC_c4;
res_c4_C_L_bsiter = c4_L - new_c4_C_L_bsiter;

v_eff_c4_C(ibs) = length(d_C_c4) - trace(CCC4_P*CC4*inv(C_std_C_c4)); % effective degrees of freedom

% Sine
d_C_s4 = [s4_L;
          Es_N_moho_7km';
          ];

% invert for MCs
std_C_s4 = [s4_std_L;
            ones(nlayersPconstr_E,1)/DAMP_E_N_moho_7km;
           ];
C_std_C_s4 = diag(std_C_s4.^2) ; % std matrix
CC4 = inv(CCC4_P'* inv(C_std_C_s4)*CCC4_P + Q_C) * CCC4_P';
MC_s4 = CC4*inv(C_std_C_s4)*d_C_s4;

new_s4_C_L = CCC4_allper*MC_s4;
res_s4_C_L = s4_L_save' - new_s4_C_L;
new_s4_C_L_bsiter = CCC4*MC_s4;
res_s4_C_L_bsiter = s4_L - new_s4_C_L_bsiter;

v_eff_s4_C(ibs) = length(d_C_s4) - trace(CCC4_P*CC4*inv(C_std_C_s4)); % effective degrees of freedom

v_eff_C(ibs) = v_eff_s4_C(ibs) + v_eff_c4_C(ibs); % effective degrees of freedom

%% Calculate chi-square
% For all data
chi2_GBH_R0(ibs) = sum(res_c2_GBH_R0.^2./c2_std_R0.^2') + sum(res_s2_GBH_R0.^2./s2_std_R0.^2'); 
chi2_GBH_R1(ibs) = sum(res_c2_GBH_R1.^2./c2_std_R1.^2') + sum(res_s2_GBH_R1.^2./s2_std_R1.^2');
chi2_G_L(ibs) = sum(res_c2_G_L.^2./c2_std_L_save.^2') + sum(res_s2_G_L.^2./s2_std_L_save.^2');
chi2_GBH(ibs) = chi2_GBH_R0(ibs) + chi2_GBH_R1(ibs) + chi2_G_L(ibs);
v_GBH(ibs) = ndata_S0*2 + ndata_S1*2 + ndata_T*2; % degrees of freedom
chi2_GBH_red(ibs) = chi2_GBH(ibs) / v_GBH(ibs);

chi2_C(ibs) = sum(res_c4_C_L.^2./c4_std_L_save.^2') + sum(res_s4_C_L.^2./s4_std_L_save.^2');
v_C(ibs) = ndata_T*2; % degrees of freedom
chi2_C_red(ibs) = chi2_C(ibs) / v_eff_C(ibs); %/ v_C(ibs);

% For bootstrap subset of data
chi2_GBH_R0_bsiter(ibs) = sum(res_c2_GBH_R0_bsiter.^2./c2_std_R0(I_bs_S0).^2') + sum(res_s2_GBH_R0_bsiter.^2./s2_std_R0(I_bs_S0).^2'); 
chi2_GBH_R1_bsiter(ibs) = sum(res_c2_GBH_R1_bsiter.^2./c2_std_R1(I_bs_S1).^2') + sum(res_s2_GBH_R1_bsiter.^2./s2_std_R1(I_bs_S1).^2');
chi2_G_L_bsiter(ibs) = sum(res_c2_G_L_bsiter.^2./c2_std_L.^2) + sum(res_s2_G_L_bsiter.^2./s2_std_L.^2);
chi2_GBH_bsiter(ibs) = chi2_GBH_R0_bsiter(ibs) + chi2_GBH_R1_bsiter(ibs) + chi2_G_L_bsiter(ibs);
v_GBH_bsiter(ibs) = length(I_bs_S0)*2 + length(I_bs_S1)*2 + length(I_bs_T2thet)*2; % degrees of freedom
chi2_GBH_red_bsiter(ibs) = chi2_GBH_bsiter(ibs) / v_GBH_bsiter(ibs);

chi2_C_bsiter(ibs) = sum(res_c4_C_L_bsiter.^2./c4_std_L.^2) + sum(res_s4_C_L_bsiter.^2./s4_std_L.^2);
v_C_bsiter(ibs) = length(I_bs_T4thet)*2; % degrees of freedom
chi2_C_red_bsiter(ibs) = chi2_C_bsiter(ibs) / v_eff_C(ibs); %/ v_C_bsiter(ibs);

%% get strength and direction from MGs and MGc =================

% G
strength_G = sqrt(MG_c2.^2+MG_s2.^2); % 2*dV/V = G/L
%strength = 2*strength ./ normalize_strength(qeqeqe:end);  
% strength_G = strength_G ;%./ L_norm;
fastdir_G = 0.5*atan2d(MG_s2,MG_c2);
fastdir_vec = [];
for idep = 1:length(depthlayer_G)
    phi_patty = 78;
    fastdir_vec(1) = fastdir_G(idep);
    fastdir_vec(2) = fastdir_G(idep)+180;
    fastdir_vec(3) = fastdir_G(idep)-180;
    [~, I] = min(abs(fastdir_vec-phi_patty));
    fastdir2_G(idep) = fastdir_vec(I);
       
end
ind = find(fastdir2_G<0);
fastdir2_G(ind) = fastdir2_G(ind)+180;

% B
strength_B = sqrt(MB_c2.^2+MB_s2.^2); % 2*dV/V = B/A
fastdir_B = 0.5*atan2d(MB_s2,MB_c2);
for idep = 1:length(depthlayer_B)
    phi_patty = 78;
    fastdir_vec(1) = fastdir_B(idep);
    fastdir_vec(2) = fastdir_B(idep)+180;
    fastdir_vec(3) = fastdir_B(idep)-180;
    [~, I] = min(abs(fastdir_vec-phi_patty));
    fastdir2_B(idep) = fastdir_vec(I);
       
end
ind = find(fastdir2_B<0);
fastdir2_B(ind) = fastdir2_B(ind)+180;

% H
strength_H = sqrt(MH_c2.^2+MH_s2.^2); % 2*dV/V = H/F
fastdir_H = 0.5*atan2d(MH_s2,MH_c2);
for idep = 1:length(depthlayer_H)
    phi_patty = 78;
    fastdir_vec(1) = fastdir_H(idep);
    fastdir_vec(2) = fastdir_H(idep)+180;
    fastdir_vec(3) = fastdir_H(idep)-180;
    [~, I] = min(abs(fastdir_vec-phi_patty));
    fastdir2_H(idep) = fastdir_vec(I);
       
end
ind = find(fastdir2_H<0);
fastdir2_H(ind) = fastdir2_H(ind)+180;
fastdir2_H(fastdir2_H<fastdir2_G) = fastdir2_H(fastdir2_H<fastdir2_G)+180; % Fix H direction


% C
strength_C = sqrt(MC_c4.^2+MC_s4.^2);
fastdir_C = (1/4)*atan2d(-MC_s4,-MC_c4);
ind = find(fastdir_C<0);
fastdir_C(ind) = fastdir_C(ind)+180;
for idep = 1:length(depthlayer_C)
    phi_patty = 78;
    fastdir_vec(1) = fastdir_C(idep);
    fastdir_vec(2) = fastdir_C(idep)+90;
    fastdir_vec(3) = fastdir_C(idep)+180;
    fastdir_vec(4) = fastdir_C(idep)+270;
    fastdir_vec(5) = fastdir_C(idep)-90;
    fastdir_vec(6) = fastdir_C(idep)-180;
    fastdir_vec(7) = fastdir_C(idep)-270;
    %[~, I] = min(abs(fastdir_vec-phi_patty+45));
    %[~, I] = min(abs(fastdir_vec-phi_patty-45));
    [~, I] = min(abs(fastdir_vec-phi_patty-45));
    fastdir2_C(idep) = fastdir_vec(I);       
end

% Calc phi and A
new_phi2_GBH_R0 = 0.5*atan2d(new_s2_GBH_R0,new_c2_GBH_R0);
new_phi2_GBH_R1 = 0.5*atan2d(new_s2_GBH_R1,new_c2_GBH_R1);
new_phi2_G_L = 0.5*atan2d(new_s2_G_L,new_c2_G_L);
new_phi4_C_L = (1/4)*atan2d(new_s4_C_L,new_c4_C_L);
ind_GBH_R0 = find(new_phi2_GBH_R0<1);
new_phi2_GBH_R0(ind_GBH_R0) = new_phi2_GBH_R0(ind_GBH_R0)+180;
ind_GBH_R1 = find(new_phi2_GBH_R1<1);
new_phi2_GBH_R1(ind_GBH_R1) = new_phi2_GBH_R1(ind_GBH_R1)+180;
ind_G_L = find(new_phi2_G_L<1);
new_phi2_G_L(ind_G_L) = new_phi2_G_L(ind_G_L)+180;
for iper = 1:length(T0periods)
    phi_patty = 78;
    phi4_vec(1) = phi4_L(iper);
    phi4_vec(2) = phi4_L(iper)+90;
    phi4_vec(3) = phi4_L(iper)+180;
    phi4_vec(4) = phi4_L(iper)+270;
    phi4_vec(5) = phi4_L(iper)-90;
    phi4_vec(6) = phi4_L(iper)-180;
    phi4_vec(7) = phi4_L(iper)-270;
   % [~, I] = min(abs(phi4_vec-phi_patty+45));
    [~, I] = min(abs(phi4_vec-phi_patty-45));
    phi4_L2(iper) = phi4_vec(I);
    
    new_phi4_vec(1) = new_phi4_C_L(iper);
    new_phi4_vec(2) = new_phi4_C_L(iper)+90;
    new_phi4_vec(3) = new_phi4_C_L(iper)+180;
    new_phi4_vec(4) = new_phi4_C_L(iper)+270;
    new_phi4_vec(5) = new_phi4_C_L(iper)-90;
    new_phi4_vec(6) = new_phi4_C_L(iper)-180;
    new_phi4_vec(7) = new_phi4_C_L(iper)-270;
%     [~, I] = min(abs(new_phi4_vec-phi_patty+45));
    [~, I] = min(abs(new_phi4_vec-phi4_L2(iper)));
    new_phi4_C_L2(iper) = new_phi4_vec(I);    
end

new_A2_GBH_R0 = sqrt(new_c2_GBH_R0.^2 + new_s2_GBH_R0.^2);
new_A2_GBH_R1 = sqrt(new_c2_GBH_R1.^2 + new_s2_GBH_R1.^2);
new_A2_G_L = sqrt(new_c2_G_L.^2 + new_s2_G_L.^2);
new_A4_C_L = sqrt(new_c4_C_L.^2 + new_s4_C_L.^2);

% Save bootstrap values
% bs_inds_S(:,ibs) = I_bs_S(:);
% bs_inds_T(:,ibs) = I_bs_T(:);
new_c2_GBH_R0_bs_all(:,ibs) = new_c2_GBH_R0(:);
new_s2_GBH_R0_bs_all(:,ibs) = new_s2_GBH_R0(:);
new_phi2_GBH_R0_bs_all(:,ibs) = new_phi2_GBH_R0(:);
new_A2_GBH_R0_bs_all(:,ibs) = new_A2_GBH_R0(:);
new_c2_GBH_R1_bs_all(:,ibs) = new_c2_GBH_R1(:);
new_s2_GBH_R1_bs_all(:,ibs) = new_s2_GBH_R1(:);
new_phi2_GBH_R1_bs_all(:,ibs) = new_phi2_GBH_R1(:);
new_A2_GBH_R1_bs_all(:,ibs) = new_A2_GBH_R1(:);
new_c2_G_L_bs_all(:,ibs) = new_c2_G_L(:);
new_s2_G_L_bs_all(:,ibs) = new_s2_G_L(:);
new_phi2_G_L_bs_all(:,ibs) = new_phi2_G_L(:);
new_A2_G_L_bs_all(:,ibs) = new_A2_G_L(:);
new_c4_C_L_bs_all(:,ibs) = new_c4_C_L(:);
new_s4_C_L_bs_all(:,ibs) = new_s4_C_L(:);
new_phi4_C_L_bs_all(:,ibs) = new_phi4_C_L2(:);
new_A4_C_L_bs_all(:,ibs) = new_A4_C_L(:);

MG_c2_bs(:,ibs) = MG_c2(:);
MB_c2_bs(:,ibs) = MB_c2(:);
MH_c2_bs(:,ibs) = MH_c2(:);
MC_c4_bs(:,ibs) = MC_c4(:);
MG_s2_bs(:,ibs) = MG_s2(:);
MB_s2_bs(:,ibs) = MB_s2(:);
MH_s2_bs(:,ibs) = MH_s2(:);
MC_s4_bs(:,ibs) = MC_s4(:);
strength_G_bs(:,ibs) = strength_G(:);
strength_B_bs(:,ibs) = strength_B(:);
strength_H_bs(:,ibs) = strength_H(:);
strength_C_bs(:,ibs) = strength_C(:);
fastdir_G_bs(:,ibs) = fastdir2_G(:);
fastdir_B_bs(:,ibs) = fastdir2_B(:);
fastdir_H_bs(:,ibs) = fastdir2_H(:);
fastdir_C_bs(:,ibs) = fastdir2_C(:);

end % end ibs

%%

new_c2_GBH_R0_bs = new_c2_GBH_R0_bs_all(:,chi2_GBH_red < chi2red_thresh);
new_s2_GBH_R0_bs = new_s2_GBH_R0_bs_all(:,chi2_GBH_red < chi2red_thresh);
new_phi2_GBH_R0_bs = new_phi2_GBH_R0_bs_all(:,chi2_GBH_red < chi2red_thresh);
new_A2_GBH_R0_bs = new_A2_GBH_R0_bs_all(:,chi2_GBH_red < chi2red_thresh);
new_c2_GBH_R1_bs = new_c2_GBH_R1_bs_all(:,chi2_GBH_red < chi2red_thresh);
new_s2_GBH_R1_bs = new_s2_GBH_R1_bs_all(:,chi2_GBH_red < chi2red_thresh);
new_phi2_GBH_R1_bs = new_phi2_GBH_R1_bs_all(:,chi2_GBH_red < chi2red_thresh);
new_A2_GBH_R1_bs = new_A2_GBH_R1_bs_all(:,chi2_GBH_red < chi2red_thresh);
new_c2_G_L_bs = new_c2_G_L_bs_all(:,chi2_GBH_red < chi2red_thresh);
new_s2_G_L_bs = new_s2_G_L_bs_all(:,chi2_GBH_red < chi2red_thresh);
new_phi2_G_L_bs = new_phi2_G_L_bs_all(:,chi2_GBH_red < chi2red_thresh);
new_A2_G_L_bs = new_A2_G_L_bs_all(:,chi2_GBH_red < chi2red_thresh);
new_c4_C_L_bs = new_c4_C_L_bs_all(:,chi2_C_red < chi2red_thresh);
new_s4_C_L_bs = new_s4_C_L_bs_all(:,chi2_C_red < chi2red_thresh);
new_phi4_C_L_bs = new_phi4_C_L_bs_all(:,chi2_C_red < chi2red_thresh);
new_A4_C_L_bs = new_A4_C_L_bs_all(:,chi2_C_red < chi2red_thresh);

new_c2_GBH_R0_mean = mean(new_c2_GBH_R0_bs,2);
new_c2_GBH_R0_std = std(new_c2_GBH_R0_bs,0,2);
new_s2_GBH_R0_mean = mean(new_s2_GBH_R0_bs,2);
new_s2_GBH_R0_std = std(new_s2_GBH_R0_bs,0,2);
new_A2_GBH_R0_mean = mean(new_A2_GBH_R0_bs,2);
new_A2_GBH_R0_med = median(new_A2_GBH_R0_bs,2);
new_A2_GBH_R0_std = std(new_A2_GBH_R0_bs,0,2);
new_phi2_GBH_R0_mean = mean(new_phi2_GBH_R0_bs,2);
new_phi2_GBH_R0_med = median(new_phi2_GBH_R0_bs,2);
new_phi2_GBH_R0_std = std(new_phi2_GBH_R0_bs,0,2);
new_c2_GBH_R1_mean = mean(new_c2_GBH_R1_bs,2);
new_c2_GBH_R1_std = std(new_c2_GBH_R1_bs,0,2);
new_s2_GBH_R1_mean = mean(new_s2_GBH_R1_bs,2);
new_s2_GBH_R1_std = std(new_s2_GBH_R1_bs,0,2);
new_A2_GBH_R1_mean = mean(new_A2_GBH_R1_bs,2);
new_A2_GBH_R1_med = median(new_A2_GBH_R1_bs,2);
new_A2_GBH_R1_std = std(new_A2_GBH_R1_bs,0,2);
new_phi2_GBH_R1_mean = mean(new_phi2_GBH_R1_bs,2);
new_phi2_GBH_R1_med = median(new_phi2_GBH_R1_bs,2);
new_phi2_GBH_R1_std = std(new_phi2_GBH_R1_bs,0,2);
new_c2_G_L_mean = mean(new_c2_G_L_bs,2);
new_c2_G_L_std = std(new_c2_G_L_bs,0,2);
new_s2_G_L_mean = mean(new_s2_G_L_bs,2);
new_s2_G_L_std = std(new_s2_G_L_bs,0,2);
new_A2_G_L_mean = mean(new_A2_G_L_bs,2);
new_A2_G_L_med = median(new_A2_G_L_bs,2);
new_A2_G_L_std = std(new_A2_G_L_bs,0,2);
new_phi2_G_L_mean = mean(new_phi2_G_L_bs,2);
new_phi2_G_L_med = median(new_phi2_G_L_bs,2);
new_phi2_G_L_std = std(new_phi2_G_L_bs,0,2);
new_c4_C_L_mean = mean(new_c4_C_L_bs,2);
new_c4_C_L_std = std(new_c4_C_L_bs,0,2);
new_s4_C_L_mean = mean(new_s4_C_L_bs,2);
new_s4_C_L_std = std(new_s4_C_L_bs,0,2);
new_A4_C_L_mean = mean(new_A4_C_L_bs,2);
new_A4_C_L_med = median(new_A4_C_L_bs,2);
new_A4_C_L_std = std(new_A4_C_L_bs,0,2);
new_phi4_C_L_mean = mean(new_phi4_C_L_bs,2);
new_phi4_C_L_med = median(new_phi4_C_L_bs,2);
new_phi4_C_L_std = std(new_phi4_C_L_bs,0,2);



MG_c2_mean = mean(MG_c2_bs,2);
MG_c2_std = std(MG_c2_bs,0,2);
MG_c2_med = median(MG_c2_bs(:,chi2_GBH_red < chi2red_thresh),2);
MG_s2_mean = mean(MG_s2_bs,2);
MG_s2_std = std(MG_s2_bs,0,2);
MG_s2_med = median(MG_s2_bs(:,chi2_GBH_red < chi2red_thresh),2);
MB_c2_mean = mean(MB_c2_bs,2);
MB_c2_std = std(MB_c2_bs,0,2);
MB_c2_med = median(MB_c2_bs(:,chi2_GBH_red < chi2red_thresh),2);
MB_s2_mean = mean(MB_s2_bs,2);
MB_s2_std = std(MB_s2_bs,0,2);
MB_s2_med = median(MB_s2_bs(:,chi2_GBH_red < chi2red_thresh),2);
MH_c2_mean = mean(MH_c2_bs,2);
MH_c2_std = std(MH_c2_bs,0,2);
MH_c2_med = median(MH_c2_bs(:,chi2_GBH_red < chi2red_thresh),2);
MH_s2_mean = mean(MH_s2_bs,2);
MH_s2_std = std(MH_s2_bs,0,2);
MH_s2_med = median(MH_s2_bs(:,chi2_GBH_red < chi2red_thresh),2);
MC_c4_mean = mean(MC_c4_bs,2);
MC_c4_std = std(MC_c4_bs,0,2);
MC_c4_med = median(MC_c4_bs(:,chi2_C_red < chi2red_thresh),2);
MC_s4_mean = mean(MC_s4_bs,2);
MC_s4_std = std(MC_s4_bs,0,2);
MC_s4_med = median(MC_s4_bs(:,chi2_C_red < chi2red_thresh),2);
% strength_bins = [0:(.08-0)/100:.07];
% fastdir_bins = [50:(150-50)/100:150];
strength_G_mean = mean(strength_G_bs,2);
strength_G_std = std(strength_G_bs,0,2);
[strength_G_cnts,strength_G_bins] = hist(strength_G_bs',nbins);
[strength_G_med, strength_G_u95, strength_G_l95, strength_G_u68, strength_G_l68 ] = get_95_68_prctile(strength_G_bs(:,chi2_GBH_red < chi2red_thresh)');
strength_B_mean = mean(strength_B_bs,2);
strength_B_std = std(strength_B_bs,0,2);
[strength_B_cnts,strength_B_bins] = hist(strength_B_bs',nbins);
[strength_B_med, strength_B_u95, strength_B_l95, strength_B_u68, strength_B_l68 ] = get_95_68_prctile(strength_B_bs(:,chi2_GBH_red < chi2red_thresh)');
strength_H_mean = mean(strength_H_bs,2);
strength_H_std = std(strength_H_bs,0,2);
[strength_H_cnts,strength_H_bins] = hist(strength_H_bs',nbins);
[strength_H_med, strength_H_u95, strength_H_l95, strength_H_u68, strength_H_l68 ] = get_95_68_prctile(strength_H_bs(:,chi2_GBH_red < chi2red_thresh)');
strength_C_mean = mean(strength_C_bs,2);
strength_C_std = std(strength_C_bs,0,2);
[strength_C_cnts,strength_C_bins] = hist(strength_C_bs',strength_G_bins);
[strength_C_med, strength_C_u95, strength_C_l95, strength_C_u68, strength_C_l68 ] = get_95_68_prctile(strength_C_bs(:,chi2_C_red < chi2red_thresh)');
fastdir_G_mean = mean(fastdir_G_bs,2);
fastdir_G_std = std(fastdir_G_bs,0,2);
[fastdir_G_med, fastdir_G_u95, fastdir_G_l95, fastdir_G_u68, fastdir_G_l68 ] = get_95_68_prctile(fastdir_G_bs(:,chi2_GBH_red < chi2red_thresh)');
[fastdir_G_cnts,fastdir_G_bins] = hist(fastdir_G_bs',nbins);
fastdir_B_mean = mean(fastdir_B_bs,2);
fastdir_B_std = std(fastdir_B_bs,0,2);
[fastdir_B_cnts,fastdir_B_bins] = hist(fastdir_B_bs',nbins);
[fastdir_B_med, fastdir_B_u95, fastdir_B_l95, fastdir_B_u68, fastdir_B_l68 ] = get_95_68_prctile(fastdir_B_bs(:,chi2_GBH_red < chi2red_thresh)');
fastdir_H_mean = mean(fastdir_H_bs,2);
fastdir_H_std = std(fastdir_H_bs,0,2);
[fastdir_H_cnts,fastdir_H_bins] = hist(fastdir_H_bs',nbins);
[fastdir_H_med, fastdir_H_u95, fastdir_H_l95, fastdir_H_u68, fastdir_H_l68 ] = get_95_68_prctile(fastdir_H_bs(:,chi2_GBH_red < chi2red_thresh)');
fastdir_C_mean = mean(fastdir_C_bs,2);
fastdir_C_std = std(fastdir_C_bs,0,2);
[fastdir_C_cnts,fastdir_C_bins] = hist(fastdir_C_bs',fastdir_G_bins);
[fastdir_C_med, fastdir_C_u95, fastdir_C_l95, fastdir_C_u68, fastdir_C_l68 ] = get_95_68_prctile(fastdir_C_bs(:,chi2_C_red < chi2red_thresh)');

new_c2_G_R0_mean = GGG2_R0_allper*MG_c2_med;
new_c2_B_R0_mean = BBB2_R0_allper*MB_c2_med;
new_c2_H_R0_mean = HHH2_R0_allper*MH_c2_med;
new_c2_GBH_R0_mean = new_c2_G_R0_mean + new_c2_B_R0_mean + new_c2_H_R0_mean;
new_c2_G_R1_mean = GGG2_R1_allper*MG_c2_med;
new_c2_B_R1_mean = BBB2_R1_allper*MB_c2_med;
new_c2_H_R1_mean = HHH2_R1_allper*MH_c2_med;
new_c2_GBH_R1_mean = new_c2_G_R1_mean + new_c2_B_R1_mean + new_c2_H_R1_mean;
new_c2_G_L_mean = GGG2_L_allper*MG_c2_med;
new_c4_C_L_mean = CCC4_allper*MC_c4_med;
new_s2_G_R0_mean = GGG2_R0_allper*MG_s2_med;
new_s2_B_R0_mean = BBB2_R0_allper*MB_s2_med;
new_s2_H_R0_mean = HHH2_R0_allper*MH_s2_med;
new_s2_GBH_R0_mean = new_s2_G_R0_mean + new_s2_B_R0_mean + new_s2_H_R0_mean;
new_s2_G_R1_mean = GGG2_R1_allper*MG_s2_med;
new_s2_B_R1_mean = BBB2_R1_allper*MB_s2_med;
new_s2_H_R1_mean = HHH2_R1_allper*MH_s2_med;
new_s2_GBH_R1_mean = new_s2_G_R1_mean + new_s2_B_R1_mean + new_s2_H_R1_mean;
new_s2_G_L_mean = GGG2_L_allper*MG_s2_med;
new_s4_C_L_mean = CCC4_allper*MC_s4_med;

% Calc phi and A
new_phi2_GBH_R0_mean = 0.5*atan2d(new_s2_GBH_R0_mean,new_c2_GBH_R0_mean);
new_phi2_GBH_R1_mean = 0.5*atan2d(new_s2_GBH_R1_mean,new_c2_GBH_R1_mean);
new_phi2_G_L_mean = 0.5*atan2d(new_s2_G_L_mean,new_c2_G_L_mean);
new_phi4_C_L_mean = (1/4)*atan2d(new_s4_C_L_mean,new_c4_C_L_mean);
ind_GBH_R0 = find(new_phi2_GBH_R0_mean<1);
new_phi2_GBH_R0_mean(ind_GBH_R0) = new_phi2_GBH_R0_mean(ind_GBH_R0)+180;
ind_GBH_R1 = find(new_phi2_GBH_R1_mean<1);
new_phi2_GBH_R1_mean(ind_GBH_R1) = new_phi2_GBH_R1_mean(ind_GBH_R1)+180;
ind_G_L = find(new_phi2_G_L_mean<1);
new_phi2_G_L_mean(ind_G_L) = new_phi2_G_L_mean(ind_G_L)+180;
phi4_vec = [];
new_phi4_vec = [];
for iper = 1:length(T0periods)
    phi_patty = 78;
    phi4_vec(1) = phi4_L(iper);
    phi4_vec(2) = phi4_L(iper)+90;
    phi4_vec(3) = phi4_L(iper)+180;
    phi4_vec(4) = phi4_L(iper)+270;
    phi4_vec(5) = phi4_L(iper)-90;
    phi4_vec(6) = phi4_L(iper)-180;
    phi4_vec(7) = phi4_L(iper)-270;
   % [~, I] = min(abs(phi4_vec-phi_patty+45));
    [~, I] = min(abs(phi4_vec-phi_patty-45));
    phi4_L2(iper) = phi4_vec(I);
    
    new_phi4_vec(1) = new_phi4_C_L_mean(iper);
    new_phi4_vec(2) = new_phi4_C_L_mean(iper)+90;
    new_phi4_vec(3) = new_phi4_C_L_mean(iper)+180;
    new_phi4_vec(4) = new_phi4_C_L_mean(iper)+270;
    new_phi4_vec(5) = new_phi4_C_L_mean(iper)-90;
    new_phi4_vec(6) = new_phi4_C_L_mean(iper)-180;
    new_phi4_vec(7) = new_phi4_C_L_mean(iper)-270;
%     [~, I] = min(abs(new_phi4_vec-phi_patty+45));
    [~, I] = min(abs(new_phi4_vec-phi4_L2(iper)));
    new_phi4_C_L2(iper) = new_phi4_vec(I);    
end
new_phi4_C_L_mean = new_phi4_C_L2;

new_A2_GBH_R0_mean = sqrt(new_c2_GBH_R0_mean.^2 + new_s2_GBH_R0_mean.^2);
new_A2_GBH_R1_mean = sqrt(new_c2_GBH_R1_mean.^2 + new_s2_GBH_R1_mean.^2);
new_A2_G_L_mean = sqrt(new_c2_G_L_mean.^2 + new_s2_G_L_mean.^2);
new_A4_C_L_mean = sqrt(new_c4_C_L_mean.^2 + new_s4_C_L_mean.^2);



%% Get R0 and R1 measurements we need to plot
I_R0 = aniso_R.mode_br == 0;
I_R1 = aniso_R.mode_br == 1;

phi2_R0 = aniso_R.phi2(I_R0);
err_phi2_R0 = err_phi2_R(I_R0);
A2_R0 = A2_R_perc(I_R0);
err_A2_R0 = err_A2_R_perc(I_R0);
c_iso_R0 = c_iso_R(I_R0);
A2_R0_perc = A2_R0;
err_A2_R0_perc = err_A2_R0;

phi2_R1 = aniso_R.phi2(I_R1);
err_phi2_R1 = aniso_R.err_phi2(I_R1);
A2_R1 = A2_R_perc(I_R1);
err_A2_R1 = err_A2_R_perc(I_R1);
c_iso_R1 = c_iso_R(I_R1);
A2_R1_perc = A2_R1;
err_A2_R1_perc = err_A2_R1;

%% Plot histograms of chi^2
f53 = figure(53); clf;
set(gcf,'position',[27           1        1185         664]);

[G_cnt, G_bins] = plot_hist([2,2,1],chi2_GBH_red,nbs);
title('G','fontsize',18);
[G_cnt_bsiter, G_bins_bsiter] = plot_hist([2,2,2],chi2_GBH_red_bsiter,nbs);
title('G (bootstraps)','fontsize',18);
[E_cnt, E_bins] = plot_hist([2,2,3],chi2_C_red,nbs);
title('E','fontsize',18);
[E_cnt_bsiter, E_bins_bsiter] = plot_hist([2,2,4],chi2_C_red_bsiter,nbs);
title('E (bootstraps)','fontsize',18);


%% Plot sine and cosine components

figure(3); clf;
% set(gcf,'position',[135   379   934   326]);
% set(gcf,'position',[135     1   934   704]);

% STRENGTH G SINE & COSINE
Nx = 2; Ny = 1;
sidegap = 0.10; topgap = 0.10; botgap = 0.10; vgap = 0.05; hgap = 0.05; cbar_bot = 0.04;
width = (1 - vgap*(Nx-1)-2*sidegap)/Nx; height = (1 - topgap - botgap - (Ny-1)*hgap)/Ny;
set(gcf,'position',[ 118   290   1000   750]);
set(gcf,'color','w');
LBLFNT=18
ix = 1; iy = 1;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height]);
hold on; set(gca,'fontsize',18);
plot(MG_c2_bs*100,depthlayer_G,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
plot(MG_s2_bs*100,depthlayer_G,'--','color',[0.7 0.85 0.7],'linewidth',2); hold on;
plot(MC_c4_bs*100,depthlayer_C,'-','color',[0.7 0.7 0.85],'linewidth',2); hold on;
plot(MC_s4_bs*100,depthlayer_C,'--','color',[0.7 0.7 0.85],'linewidth',2); hold on;
h1(2) = plot(MB_c2_mean*100,depthlayer_B,'-','color',[1 0.7 0.7],'linewidth',4);hold on;
h1(3) = plot(MH_c2_mean*100,depthlayer_H,'-','color',[0.7 0.7 1],'linewidth',4);hold on;
h1(1) = plot(MG_c2_mean*100,depthlayer_G,'-','color',[0 0.7 0],'linewidth',4);hold on;
h1(4) = plot(MC_c4_mean*100,depthlayer_C,'-','color',[0 0.4470 0.7410],'linewidth',4);hold on;
plot(MB_s2_mean*100, depthlayer_B,'--','color',[1 0.7 0.7],'linewidth',4);hold on;
plot(MH_s2_mean*100, depthlayer_H,'--','color',[0.7 0.7 1],'linewidth',4);hold on;
plot(MG_s2_mean*100, depthlayer_G,'--','color',[0 0.7 0],'linewidth',4);hold on;
plot(MC_s4_mean*100, depthlayer_C,'--','color',[0 0.4470 0.7410],'linewidth',4);hold on;
% xlim([0 10]); 
% ylim(ylims);

ylabel('Period (s)','fontsize',18)
xlabel('Strength (%)','fontsize',18)
set(gca,'box','on','LineWidth', 1.5,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
legend(h1,{'G/L (2\theta)','B/A (2\theta)','H/F (2\theta)','E/N (4\theta)'},'fontsize',18,'location','southwest','box','off');
title('Model:  -Cos  --Sin');

set(gca,'YDir','reverse')

% STRENGTH DATA SINE & COSINE
ix = 1;iy = 2;
left = sidegap + (iy-1)*(vgap+width);bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height]);
set(gca,'fontsize',18);

semilogy(new_c2_GBH_R0_bs*2*100, S0periods,'-','color',[0.7 0.85 0.7],'linewidth',2);hold on;
plot(new_c2_GBH_R1_bs*2*100, S1periods,'-','color',[0.7 0.85 0.7],'linewidth',2);hold on;
plot(new_c2_G_L_bs*2*100, T0periods,'-','color',[0.7 0.85 0.7],'linewidth',2);hold on;
plot(new_c4_C_L_bs*2*100, T0periods,'-','color',[0.7 0.7 0.85],'linewidth',2);hold on;
plot(new_s2_GBH_R0_bs*2*100, S0periods,'--','color',[0.7 0.85 0.7],'linewidth',2);hold on;
plot(new_s2_GBH_R1_bs*2*100, S1periods,'--','color',[0.7 0.85 0.7],'linewidth',2);hold on;
plot(new_s2_G_L_bs*2*100, T0periods,'--','color',[0.7 0.85 0.7],'linewidth',2);hold on;
plot(new_s4_C_L_bs*2*100, T0periods,'--','color',[0.7 0.7 0.85],'linewidth',2);hold on;

plot(new_c2_GBH_R0_mean*2*100, S0periods,'-','color',[0 0.7 0],'linewidth',4);hold on;
h = errorbar(c2_R0*2*100,S0periods,c2_std_R0*2*100,'horizontal','x','markersize',12,'linewidth',2,'color','k');
h.CapSize = 0;
plot(new_s2_GBH_R0_mean*2*100, S0periods,'--','color',[0 0.7 0],'linewidth',4);hold on;
h = errorbar(s2_R0*2*100,S0periods,s2_std_R0*2*100,'horizontal','x','markersize',12,'linewidth',2,'color','k');
h.CapSize = 0;

plot(new_c2_GBH_R1_mean*2*100, S1periods,'-','color',[0 0.7 0],'linewidth',4);hold on;
h = errorbar(c2_R1*2*100,S1periods,c2_std_R1*2*100,'horizontal','x','markersize',12,'linewidth',2,'color','k');
h.CapSize = 0;
plot(new_s2_GBH_R1_mean*2*100, S1periods,'--','color',[0 0.7 0],'linewidth',4);hold on;
h = errorbar(s2_R1*2*100,S1periods,s2_std_R1*2*100,'horizontal','x','markersize',12,'linewidth',2,'color','k');
h.CapSize = 0;

plot(new_c2_G_L_mean*2*100, T0periods,'-','color',[0 0.4 0],'linewidth',4);hold on;
h = errorbar(c2_L_save*2*100,T0periods,c2_std_L_save*2*100,'horizontal','x','markersize',12,'linewidth',2,'color','k');
h.CapSize = 0;
plot(new_s2_G_L_mean*2*100, T0periods,'--','color',[0 0.4 0],'linewidth',4);hold on;
h = errorbar(s2_L_save*2*100,T0periods,s2_std_L_save*2*100,'horizontal','x','markersize',12,'linewidth',2,'color','k');
h.CapSize = 0;

plot(new_c4_C_L_mean*2*100, T0periods,'-','color',[0 0 0.7],'linewidth',4);hold on;
h = errorbar(c4_L_save*2*100,T0periods,c4_std_L_save*2*100,'horizontal','x','markersize',12,'linewidth',2,'color','k');
h.CapSize = 0;
plot(new_s4_C_L_mean*2*100, T0periods,'--','color',[0 0 0.7],'linewidth',4);hold on;
h = errorbar(s4_L_save*2*100,T0periods,s4_std_L_save*2*100,'horizontal','x','markersize',12,'linewidth',2,'color','k');
h.CapSize = 0;
% plot([ 78 78],[0 250],'k--','Linewidth',3);
% APM = 296.74;
% plot([ APM-180 APM-180 ],[0 250],'--','color',[0.6 0.6 0.6], 'Linewidth',3);
set(gca,'YDir','reverse')
% xlim([0 180]);
ylim([min(S1periods)-1 max(S0periods)+10]);
xlabel('Strength (%)','fontsize',18)
set(gca,'box','on','LineWidth', 1.5,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% legend(h1,{'G/L (2\theta)','B/A (2\theta)','H/F (2\theta)'},'fontsize',18);
title('Data:  -Cos  --Sin');

if isfig
%     save2pdf([figpath,'azi_depth_GBHE_bootstrap_alphF',num2str(alphF_GBH,'%2.1e'),'_nperm',num2str(nperm_S),'_nbs',num2str(nbs),'_ORALS_long_DAMP0_NORM_perc_MN86_GL_sinecos_errprop_patch.pdf'],3,1000)
end


%% plot inversion result for P2P strength and fastdir as function of Depth (SPLIT) WITH VELOCITY
figure(105);clf;
shift = 0.25;
% VSV
Nx = 2.5; Ny = 1;
sidegap = 0.10; topgap = 0.10; botgap = 0.10; vgap = 0.05; hgap = 0.05; cbar_bot = 0.04;
width = (1 - vgap*(Nx-1)-2*sidegap)/Nx; height = (1 - topgap - botgap - (Ny-1)*hgap)/Ny;
%set(gcf,'position',[ 118   290   1000   750]);
set(gcf,'position',[1           1        1280         704]);
set(gcf,'color','w');
LBLFNT=18
ix = 1; iy = 1;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width*0.75,height*.75]);
hold on; set(gca,'fontsize',18,'linewidth',2);
plot(Vsv/1000,depthlayer_G,'-','color',[0 0 0],'linewidth',4);hold on;
xlim([4.2 4.8]); 
ylim([40 300]);
ylabel('Depth (km)','fontsize',18)
xlabel('Velocity (km/s)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'YDir','reverse')

ix = 1.72; iy = 1;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
axes('position',[left,bot,width*0.75,height*.3]);
hold on; set(gca,'fontsize',18,'linewidth',2);
plot(Vsv/1000,depthlayer_G,'-','color',[0 0 0],'linewidth',4);hold on;
xlim([4.2 4.8]); 
ylim([11 40]);
% ylabel('Depth (km)','fontsize',18)
% xlabel('Strength (%)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'YDir','reverse','xticklabel',[])

% STRENGTH
ix = 1; iy = 1.5+shift;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height*.75]);
hold on; set(gca,'fontsize',18,'linewidth',2);
% plot(strength_G_bs*100,depthlayer_G,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
h = patch([strength_G_l95*100 fliplr(strength_G_u95*100)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([strength_G_l68*100 fliplr(strength_G_u68*100)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
% plot(strength_C_bs*100,depthlayer_C,'-','color',[0.7 0.7 0.85],'linewidth',2); hold on;
h = patch([strength_C_l95*100 fliplr(strength_C_u95*100)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([strength_C_l68*100 fliplr(strength_C_u68*100)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
% plot(strength_B_mean*100,depthlayer_B,'-','color',[1 0.7 0.7],'linewidth',4);hold on;
% plot(strength_H_mean*100,depthlayer_H,'-','color',[0.7 0.7 1],'linewidth',4);hold on;
plot(strength_G_med*100,depthlayer_G,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);hold on;
plot(strength_C_med*100,depthlayer_C,'-','color',[0 0.4470 0.7410],'linewidth',4);hold on;
xlim([0 8]); 
ylim([40 300]);
% ylabel('Depth (km)','fontsize',18)
xlabel('Strength (%)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'YDir','reverse')

ix = 1.72; iy = 1.5+shift;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
axes('position',[left,bot,width,height*.3]);
hold on; set(gca,'fontsize',18,'linewidth',2);
% plot(strength_G_bs*100,depthlayer_G,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
h = patch([strength_G_l95*100 fliplr(strength_G_u95*100)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([strength_G_l68*100 fliplr(strength_G_u68*100)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
% plot(strength_C_bs*100,depthlayer_C,'-','color',[0.7 0.7 0.85],'linewidth',2); hold on;
h = patch([strength_C_l95*100 fliplr(strength_C_u95*100)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([strength_C_l68*100 fliplr(strength_C_u68*100)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
% plot(strength_B_mean*100,depthlayer_B,'-','color',[1 0.7 0.7],'linewidth',4);hold on;
% plot(strength_H_mean*100,depthlayer_H,'-','color',[0.7 0.7 1],'linewidth',4);hold on;
plot(strength_G_med*100,depthlayer_G,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);hold on;
plot(strength_C_med*100,depthlayer_C,'-','color',[0 0.4470 0.7410],'linewidth',4);hold on;
clr = lines(length(cij_ind));
iter = 0;
bar_len = [15, 19]; %[15, 17, 19, 19];
for ic = cij_ind
    iter = iter+1;
    strength_g = cij_calc(ic).strength_g;
    strength_c = cij_calc(ic).strength_c;
%     h(ic) = plot([strength_g strength_g],[0 20],'-','color',clr(iter,:),'linewidth',4); hold on;
%     plot([strength_c strength_c],[0 20],'--','color',clr(iter,:),'linewidth',4);
    h(ic) = plot([strength_g strength_g],[0 bar_len(iter)],'-','color',[0.2660    0.6740    0.1880],'linewidth',4); hold on;
    plot([strength_c strength_c],[0 bar_len(iter)],'-','color',[0 0.4470 0.7410],'linewidth',4);
end
xlim([0 8]); 
ylim([11 40]);
% ylabel('Depth (km)','fontsize',18)
% xlabel('Strength (%)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'YDir','reverse','xticklabel',[])


% AZIMUTH
ix = 1;iy = 2.5+shift;
left = sidegap + (iy-1)*(vgap+width);bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height*0.75]);
hold on;set(gca,'fontsize',18,'linewidth',2);hold on;set(gca,'fontsize',18,'linewidth',2);
% plot(fastdir_G_bs,depthlayer_G,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
h = patch([fastdir_G_l95 fliplr(fastdir_G_u95)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([fastdir_G_l68 fliplr(fastdir_G_u68)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
% plot(strength_C_bs*100,depthlayer_C,'-','color',[0.7 0.7 0.85],'linewidth',2); hold on;
h = patch([fastdir_C_l95 fliplr(fastdir_C_u95)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([fastdir_C_l68 fliplr(fastdir_C_u68)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
% h1(2) = plot(fastdir_B_mean, depthlayer_B,'-','color',[1 0.7 0.7],'linewidth',4);hold on;
% h1(3) = plot(fastdir_H_mean, depthlayer_H,'-','color',[0.7 0.7 1],'linewidth',4);hold on;
h1(1) = plot(fastdir_G_mean, depthlayer_G,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);hold on;
h1(2) = plot(fastdir_C_mean, depthlayer_C,'-','color',[0 0.4470 0.7410],'linewidth',4);hold on;
plot([ 78 78],[0 250],'k--','Linewidth',3);
APM = 296.74;
plot([ APM-180 APM-180 ],[0 400],'--','color',[0.6 0.6 0.6], 'Linewidth',3);
set(gca,'YDir','reverse')
xlim([30 200]);
ylim([40 300]);
xlabel('Azimuth (\circ)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
set(gca,'xtick',[0:30:180]);
% legend(h1,{'G/L (2\theta)','B/A (2\theta)','H/F (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');
% legend(h1,{'G/L (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');

ix = 1.72;iy = 2.5+shift;
left = sidegap + (iy-1)*(vgap+width);bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height*0.3]);
hold on;set(gca,'fontsize',18,'linewidth',2);hold on;set(gca,'fontsize',18,'linewidth',2);
% plot(fastdir_G_bs,depthlayer_G,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
h = patch([fastdir_G_l95 fliplr(fastdir_G_u95)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([fastdir_G_l68 fliplr(fastdir_G_u68)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
% plot(strength_C_bs*100,depthlayer_C,'-','color',[0.7 0.7 0.85],'linewidth',2); hold on;
h = patch([fastdir_C_l95 fliplr(fastdir_C_u95)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([fastdir_C_l68 fliplr(fastdir_C_u68)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
% h1(2) = plot(fastdir_B_mean, depthlayer_B,'-','color',[1 0.7 0.7],'linewidth',4);hold on;
% h1(3) = plot(fastdir_H_mean, depthlayer_H,'-','color',[0.7 0.7 1],'linewidth',4);hold on;
h1(1) = plot(fastdir_G_mean, depthlayer_G,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);hold on;
h1(2) = plot(fastdir_C_mean, depthlayer_C,'-','color',[0 0.4470 0.7410],'linewidth',4);hold on;
plot([ 78 78],[0 250],'k--','Linewidth',3);
APM = 296.74;
plot([ APM-180 APM-180 ],[0 400],'--','color',[0.6 0.6 0.6], 'Linewidth',3);
set(gca,'YDir','reverse')
xlim([30 200]);
ylim([11 40]);
% xlabel('Azimuth (\circ)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
set(gca,'xticklabel',[]);
set(gca,'xtick',[0:30:180]);
% legend(h1,{'G/L (2\theta)','B/A (2\theta)','H/F (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');
legend(h1,{'G/L (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');


% SUBPLOT AZIMUTH
ix = 1.0;
iy = 2.5+shift;

left = width/1.3+ (iy-1)*(vgap+width);
bot = 0.08 + (ix-1)*(hgap+height);
axes('position',[left+width/7,bot+height/10,width/3,height/3]);
box on
% new_phi2_GBH_R0 = 0.5*atan2d(new_s2_GBH_R0,new_c2_GBH_R0);
% new_phi2_GBH_R1 = 0.5*atan2d(new_s2_GBH_R1,new_c2_GBH_R1);
% new_phi2_G_L = 0.5*atan2d(new_s2_G_L,new_c2_G_L);
% new_phi4_C_L = (1/4)*atan2d(new_s4_C_L,new_c4_C_L);
% ind_G_R0 = find(new_phi2_GBH_R0<1);
% new_phi2_GBH_R0(ind_G_R0) = new_phi2_GBH_R0(ind_G_R0)+180;
% ind_G_R1 = find(new_phi2_GBH_R1<1);
% new_phi2_GBH_R1(ind_G_R1) = new_phi2_GBH_R1(ind_G_R1)+180;
% ind_G_L = find(new_phi2_G_L<1);
% new_phi2_G_L(ind_G_L) = new_phi2_G_L(ind_G_L)+180;
for iper = 1:length(T0periods)
    phi_patty = 78;
    phi4_vec(1) = phi4_L(iper);
    phi4_vec(2) = phi4_L(iper)+90;
    phi4_vec(3) = phi4_L(iper)+180;
    phi4_vec(4) = phi4_L(iper)+270;
    phi4_vec(5) = phi4_L(iper)-90;
    phi4_vec(6) = phi4_L(iper)-180;
    phi4_vec(7) = phi4_L(iper)-270;
   % [~, I] = min(abs(phi4_vec-phi_patty+45));
    [~, I] = min(abs(phi4_vec-phi_patty-45));
    phi4_L2(iper) = phi4_vec(I);  
end
plot(new_phi2_GBH_R0_bs,S0periods,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
plot(new_phi2_GBH_R0_mean,S0periods,'-','color',[0.2660    0.6740    0.1880],'linewidth',2)
plot(phi2_R0,S0periods,'o','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',8,'markerfacecolor',[0.2660    0.6740    0.1880]);hold on;
h = errorbar(phi2_R0,S0periods,err_phi2_R0,'horizontal','.','linewidth',1.5,'color','k');
h.CapSize = 0;
xlim([50 200]);
ylim([min(S0periods)-10 max(S0periods)+10])
plot([ 78 78],[0 250],'k--','Linewidth',1.5);
APM = 296.74;
plot([ APM-180 APM-180 ],[0 400],'--','color',[0.6 0.6 0.6], 'Linewidth',1.5);
% ylabel('Period (s)','fontsize',12);
xlabel('\psi (\circ)','fontsize',12)
set(gca,'xtick',[0 45 90 135 180],'fontsize',12,'linewidth',1.5,'ydir','reverse','TickLength',[0.03, 0.01]);

ix = 1.32;
iy = 2.5+shift;
left = width/1.3+ (iy-1)*(vgap+width);
bot = 0.08 + (ix-1)*(hgap+height);
axes('position',[left+width/7,bot+height/10,width/3,height/5]);
box on
plot(new_phi2_GBH_R1_bs,S1periods,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
plot(new_phi2_G_L_bs,T0periods,'-','color',[0.7 0.85 0.7],'linewidth',2)
plot(new_phi4_C_L_bs,T0periods,'-','color',[0.6824 0.8353 0.9765],'linewidth',2)
plot(new_phi2_GBH_R1_mean,S1periods,'--','color',[0.2660    0.6740    0.1880],'linewidth',2);
plot(new_phi2_G_L_mean,T0periods,'-','color',[0.2660    0.6740    0.1880],'linewidth',2);
plot(new_phi4_C_L_mean,T0periods,'-','color',[0 0.4470 0.7410],'linewidth',2);
plot(phi2_R1,S1periods,'o','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',8,'markerfacecolor',[0.2660    0.6740    0.1880]);hold on;
plot(phi2_L,T0periods,'^','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',8,'markerfacecolor',[0.2660    0.6740    0.1880]);hold on;
plot(phi4_L2,T0periods,'^','color',[0 0.4470 0.7410],'linewidth',2,'markersize',8,'markerfacecolor',[0 0.4470 0.7410]);hold on;
h = errorbar(phi2_R1,S1periods,err_phi2_R1,'horizontal','.','linewidth',1.5,'color','k');
h.CapSize = 0;
h = errorbar(phi2_L,T0periods,err_phi2_L,'horizontal','.','linewidth',1.5,'color','k');
h.CapSize = 0;
h = errorbar(phi4_L2,T0periods,err_phi4_L,'horizontal','.','linewidth',1.5,'color','k');
h.CapSize = 0;
plot([ 78 78],[0 250],'k--','Linewidth',1.5);
plot([ 78 78]+45,[0 250],'k--','Linewidth',1.5);
plot([ 78 78]+90,[0 250],'k--','Linewidth',1.5);
xlim([50 200]);
ylim([min(S1periods)-0.25 max(S1periods)+0.25])
ylabel('Period (s)','fontsize',12);
% xlabel('\psi (\circ)','fontsize',12)
set(gca,'xtick',[0 45 90 135 180],'xticklabel',[],'fontsize',12,'linewidth',1.5,'ydir','reverse','TickLength',[0.03, 0.01]);


% SUBPLOT STRENGTH
ix = 1;
iy = 1.5+shift;
left = width/1.3+ (iy-1)*(vgap+width);
bot = 0.08 + (ix-1)*(hgap+height);
axes('position',[left+width/7,bot+height/10,width/3,height/3]);
box on;
plot(new_A2_GBH_R0_bs*2*100,S0periods,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
plot(new_A2_GBH_R0_mean*2*100,S0periods,'-','color',[0.2660    0.6740    0.1880],'linewidth',2);
plot(A2_R0_perc*2*100,S0periods,'o','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',8,'markerfacecolor',[0.2660    0.6740    0.1880]); hold on;
h = errorbar(A2_R0_perc*2*100,S0periods,err_A2_R0_perc*2*100,'horizontal','.','linewidth',1.5,'color','k');
h.CapSize = 0;
ylim([min(S0periods)-10 max(S0periods)+10])
% ylabel('Period (s)','fontsize',12);
xlabel('2A (%)','fontsize',12);
set(gca,'fontsize',12,'linewidth',1.5,'ydir','reverse','TickLength',[0.03, 0.01])
xlim([0 5])

ix = 1.32;
iy = 1.5+shift;
left = width/1.3+ (iy-1)*(vgap+width);
bot = 0.08 + (ix-1)*(hgap+height);
axes('position',[left+width/7,bot+height/10,width/3,height/5]);
box on;
plot(new_A2_GBH_R1_bs*2*100,S1periods-.1,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
plot(new_A2_G_L_bs*2*100,T0periods+.1,'-','color',[0.7 0.85 0.7],'linewidth',2)
plot(new_A4_C_L_bs*2*100,T0periods,'-','color',[0.6824 0.8353 0.9765],'linewidth',2)
plot(new_A2_GBH_R1_mean*2*100,S1periods-.1,'--','color',[0.2660    0.6740    0.1880],'linewidth',2);
plot(new_A2_G_L_mean*2*100,T0periods+.1,'-','color',[0.2660    0.6740    0.1880],'linewidth',2);
plot(new_A4_C_L_mean*2*100,T0periods,'-','color',[0 0.4470 0.7410],'linewidth',2);
plot(A2_R1_perc*2*100,S1periods-.1,'o','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',8,'markerfacecolor',[0.2660    0.6740    0.1880]); hold on;
plot(A2_L_perc*2*100,T0periods+.1,'^','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',8,'markerfacecolor',[0.2660    0.6740    0.1880]);
plot(A4_L_perc*2*100,T0periods,'^','color',[0 0.4470 0.7410],'linewidth',2,'markersize',8,'markerfacecolor',[0 0.4470 0.7410]);
h = errorbar(A2_R1_perc*2*100,S1periods-.1,err_A2_R1_perc*2*100,'horizontal','.','linewidth',1.5,'color','k');
h.CapSize = 0;
h = errorbar(A2_L_perc*2*100,T0periods+.1,err_A2_L_perc*2*100,'horizontal','.','linewidth',1.5,'color','k');
h.CapSize = 0;
h = errorbar(A4_L_perc*2*100,T0periods,err_A4_L_perc*2*100,'horizontal','.','linewidth',1.5,'color','k');
h.CapSize = 0;
ylim([min(S1periods)-.25 max(S1periods)+.25])
ylabel('Period (s)','fontsize',12);
% xlabel('2A (%)','fontsize',12);
set(gca,'xticklabel',[],'fontsize',12,'linewidth',1.5,'ydir','reverse','TickLength',[0.03, 0.01])
xlim([0 5])

if isfig
%     save2pdf([figpath,'azi_depth_GBHE_bootstrap_alphF',num2str(alphF_GBH,'%2.1e'),'_nperm',num2str(nperm_S),'_nbs',num2str(nbs),'_ES172_long_DAMP0_NORM_perc_MN86_GL_SPLIT_errprop_patch.pdf'],103,1000)
    export_fig([figpath,'azi_depth_GBHE_bootstrap_alphF',num2str(alphF_GBH,'%2.1e'),'_nbs',num2str(nbs),'_ORALS_long_DAMP0_NORM_perc_MN86_GL_SPLIT_errprop_patch2_vel.pdf'],'-pdf','-q100','-p0.02','-painters',105)
end

%% Heat Maps 1
figure(51); clf;
set(gcf,'position',[1           1        1250         697]);
max_cbar = 0.95;
ncont = 10;
cmapp = parula;

%G
subplot(1,4,1); box on;
contourf(strength_G_bins*100,depthlayer_G,strength_G_cnts'./max(strength_G_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(strength_G_l95*100,depthlayer_G,'-y','linewidth',3);
plot(strength_G_u95*100,depthlayer_G,'-y','linewidth',3);
plot(strength_G_l68*100,depthlayer_G,'-r','linewidth',3);
plot(strength_G_u68*100,depthlayer_G,'-r','linewidth',3);
ylim([12 300]);
caxis([0 max_cbar]);
ylabel('Depth (km)','fontsize',18)
xlabel('Strength (%)','fontsize',18)
title('G/L','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');
colormap(cmapp);

subplot(1,4,2); box on;
contourf(fastdir_G_bins,depthlayer_G,fastdir_G_cnts'./max(fastdir_G_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(fastdir_G_l95,depthlayer_G,'-y','linewidth',3);
plot(fastdir_G_u95,depthlayer_G,'-y','linewidth',3);
plot(fastdir_G_l68,depthlayer_G,'-r','linewidth',3);
plot(fastdir_G_u68,depthlayer_G,'-r','linewidth',3);
% ylabel('Depth (km)','fontsize',18)
xlabel('Azimuth (\circ)','fontsize',18)
title('G/L','fontsize',18);
ylim([12 300]);
xlim([50 150]);
caxis([0 max_cbar]);
plot([ 78 78],[0 250],'--','color',[0.8 0.8 0.8],'Linewidth',3);
APM = 296.74;
plot([ APM-180 APM-180 ],[0 400],'--','color',[1 1 1], 'Linewidth',3);
set(gca,'YDir','reverse')
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');
set(gca,'xtick',[0:30:180]);

% C
subplot(1,4,3); box on;
contourf(strength_C_bins*100,depthlayer_C,strength_C_cnts'./max(strength_C_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(strength_C_l95*100,depthlayer_C,'-y','linewidth',2);
plot(strength_C_u95*100,depthlayer_C,'-y','linewidth',2);
plot(strength_C_l68*100,depthlayer_C,'-r','linewidth',2);
plot(strength_C_u68*100,depthlayer_C,'-r','linewidth',2);
% ylim([12 40]);
caxis([0 max_cbar]);
ylim([12 300])
xlim([min(strength_G_bins*100), max(strength_G_bins*100)]);
xlabel('Strength (%)','fontsize',18)
title('E/N','fontsize',18);
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');

subplot(1,4,4); box on;
contourf(fastdir_C_bins,depthlayer_C,fastdir_C_cnts'./max(fastdir_C_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(fastdir_C_l95,depthlayer_C,'-y','linewidth',2);
plot(fastdir_C_u95,depthlayer_C,'-y','linewidth',2);
plot(fastdir_C_l68,depthlayer_C,'-r','linewidth',2);
plot(fastdir_C_u68,depthlayer_C,'-r','linewidth',2);
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

%% Heat Maps 2
max_cbar = 0.95;
ncont = 50;
cmapp = parula;

figure(52); clf;
shift = 0.7;
Nx = 2.85; Ny = 1;
sidegap = 0.05; topgap = 0.10; botgap = 0.10; vgap = 0.05; hgap = 0.05; cbar_bot = 0.04;
width = (1 - vgap*(Nx-1)-2*sidegap)/Nx; height = (1 - topgap - botgap - (Ny-1)*hgap)/Ny;
%set(gcf,'position',[ 118   290   1000   750]);
set(gcf,'position',[1           1        1280         704]);
LBLFNT=18;

%G Strength
ix = 1; iy = 1;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width*0.75,height*.75]);
hold on; set(gca,'fontsize',18,'linewidth',2);
contourf(strength_G_bins*100,depthlayer_G,strength_G_cnts'./max(strength_G_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(strength_G_l95*100,depthlayer_G,'-y','linewidth',3);
plot(strength_G_u95*100,depthlayer_G,'-y','linewidth',3);
plot(strength_G_l68*100,depthlayer_G,'-r','linewidth',3);
plot(strength_G_u68*100,depthlayer_G,'-r','linewidth',3);
ylim([max(depthlayer_C) 300]);
caxis([0 max_cbar]);
ylabel('Depth (km)','fontsize',18)
xlabel('Strength (%)','fontsize',18)
title('G/L','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');
colormap(cmapp);

ix = 1.72; iy = 1;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
axes('position',[left,bot,width*0.75,height*.3]);
hold on; set(gca,'fontsize',18,'linewidth',2);
contourf(strength_G_bins*100,depthlayer_G,strength_G_cnts'./max(strength_G_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(strength_G_l95*100,depthlayer_G,'-y','linewidth',3);
plot(strength_G_u95*100,depthlayer_G,'-y','linewidth',3);
plot(strength_G_l68*100,depthlayer_G,'-r','linewidth',3);
plot(strength_G_u68*100,depthlayer_G,'-r','linewidth',3);
ylim([min(depthlayer_G) max(depthlayer_C)]);
caxis([0 max_cbar]);
ylabel('Depth (km)','fontsize',18)
title('G/L','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');
set(gca,'xticklabel',[]);
colormap(cmapp);

% G fast direction
ix = 1; iy = 1+shift;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width*0.75,height*.75]);
hold on; set(gca,'fontsize',18,'linewidth',2);
contourf(fastdir_G_bins,depthlayer_G,fastdir_G_cnts'./max(fastdir_G_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(fastdir_G_l95,depthlayer_G,'-y','linewidth',3);
plot(fastdir_G_u95,depthlayer_G,'-y','linewidth',3);
plot(fastdir_G_l68,depthlayer_G,'-r','linewidth',3);
plot(fastdir_G_u68,depthlayer_G,'-r','linewidth',3);
% ylabel('Depth (km)','fontsize',18)
xlabel('Azimuth (\circ)','fontsize',18)
title('G/L','fontsize',18);
ylim([max(depthlayer_C) 300]);
xlim([50 150]);
caxis([0 max_cbar]);
plot([ 78 78],[0 250],'--','color',[0.8 0.8 0.8],'Linewidth',3);
APM = 296.74;
plot([ APM-180 APM-180 ],[0 400],'--','color',[1 1 1], 'Linewidth',3);
set(gca,'YDir','reverse')
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');
set(gca,'xtick',[0:30:180],'yticklabel',[]);

ix = 1.72; iy = 1+shift;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
axes('position',[left,bot,width*0.75,height*.3]);
hold on; set(gca,'fontsize',18,'linewidth',2);
contourf(fastdir_G_bins,depthlayer_G,fastdir_G_cnts'./max(fastdir_G_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(fastdir_G_l95,depthlayer_G,'-y','linewidth',3);
plot(fastdir_G_u95,depthlayer_G,'-y','linewidth',3);
plot(fastdir_G_l68,depthlayer_G,'-r','linewidth',3);
plot(fastdir_G_u68,depthlayer_G,'-r','linewidth',3);
% ylabel('Depth (km)','fontsize',18)
title('G/L','fontsize',18);
ylim([min(depthlayer_G) max(depthlayer_C)]);
xlim([50 150]);
caxis([0 max_cbar]);
plot([ 78 78],[0 250],'--','color',[0.8 0.8 0.8],'Linewidth',3);
APM = 296.74;
plot([ APM-180 APM-180 ],[0 400],'--','color',[1 1 1], 'Linewidth',3);
set(gca,'YDir','reverse')
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');
set(gca,'xticklabel',[],'yticklabel',[]);

% C Strength
ix = 1.72;iy = 2+shift-0.3;
left = sidegap + (iy-1)*(vgap+width);bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width*0.75,height*0.3]);
hold on;set(gca,'fontsize',18,'linewidth',2);
contourf(strength_C_bins*100,depthlayer_C,strength_C_cnts'./max(strength_C_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(strength_C_l95*100,depthlayer_C,'-y','linewidth',2);
plot(strength_C_u95*100,depthlayer_C,'-y','linewidth',2);
plot(strength_C_l68*100,depthlayer_C,'-r','linewidth',2);
plot(strength_C_u68*100,depthlayer_C,'-r','linewidth',2);
% ylim([12 40]);
caxis([0 max_cbar]);
ylim([min(depthlayer_C) max(depthlayer_C)])
xlim([min(strength_G_bins*100), max(strength_G_bins*100)]);
xlabel('Strength (%)','fontsize',18)
title('E/N','fontsize',18);
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');
set(gca,'yticklabel',[]);

% C fast direction
ix = 1.72;iy = 3+shift-0.6;
left = sidegap + (iy-1)*(vgap+width);bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width*0.75,height*0.3]);
hold on;set(gca,'fontsize',18,'linewidth',2);
contourf(fastdir_C_bins,depthlayer_C,fastdir_C_cnts'./max(fastdir_C_cnts)',ncont,'linestyle','none'); hold on; shading flat;
plot(fastdir_C_l95,depthlayer_C,'-y','linewidth',2);
plot(fastdir_C_u95,depthlayer_C,'-y','linewidth',2);
plot(fastdir_C_l68,depthlayer_C,'-r','linewidth',2);
plot(fastdir_C_u68,depthlayer_C,'-r','linewidth',2);
% ylabel('Depth (km)','fontsize',18)
xlabel('Azimuth (\circ)','fontsize',18)
title('E/N','fontsize',18);
caxis([0 max_cbar]);
% ylim([12 40]);
ylim([min(depthlayer_C) max(depthlayer_C)])
xlim([50 150]);
plot([ 78 78],[0 250],'--','color',[0.8 0.8 0.8],'Linewidth',3);
APM = 296.74;
plot([ APM-180 APM-180 ],[0 400],'--','color',[1 1 1], 'Linewidth',3);
set(gca,'YDir','reverse')
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT,'ydir','reverse');
set(gca,'xtick',[0:30:180],'yticklabel',[]);

%% CALC CIJs
% % Gc = MG_c2_mean.*L_norm;
% % Bc = MB_c2_mean.*A_norm;
% % Hc = MH_c2_mean.*F_norm;
% % Ec = MC_c4_mean.*N_norm;
% Gc = sqrt(MG_c2_mean.^2+MG_s2_mean.^2).*L_norm/2;
% Bc = sqrt(MB_c2_mean.^2+MB_s2_mean.^2).*A_norm/2;
% Hc = sqrt(MH_c2_mean.^2+MH_s2_mean.^2).*F_norm/2;
% Ec = sqrt(MC_c4_mean.^2+MC_s4_mean.^2).*N_norm/2;
% 
% % Compare cij
% [cij_model, dep1] = ACFLNGBHE_to_cij(A_norm,C_norm,F_norm,L_norm,N_norm,Gc,Bc,Hc,Ec,depthlayer_G,20);
% cij_model = cij_model / 1e9;
% figure(5); clf;
% set(gcf,'position',[2         166        1190         539]);
% for icij = 1:length(cij)
%     subplot(2,3,icij);
%     cij_resid = (cij_model-cij(icij).c) ./ cij_model;
%     cij_resid(isnan(cij_resid)) = 0;
%     imagesc(cij_resid); hold on;
%     % plot(6,6,'ok','markerfacecolor',[0 0 0]);
%     fs = 10;
%     text(1,1,'PH','fontsize',fs);
%     text(2,2,'PH','fontsize',fs);
%     text(3,3,'PV','fontsize',fs);
%     text(4,4,'SV','fontsize',fs);
%     text(5,5,'SV','fontsize',fs);
%     text(6,6,'SH','fontsize',fs);
%     text(1.7,1,'PH,SH','fontsize',fs);
%     text(2.6,1,'\eta,PH,SV','fontsize',fs);
%     text(2.6,2,'\eta,PH,SV','fontsize',fs);
%     colorbar;
%     colormap(redblue);
%     caxis([-0.25 0.25]);
%     title(cij(icij).ref);
%     set(gca,'fontsize',16);
% end
% 
% % Compare cij misfit
% figure(6); clf;
% clr = lines(length(cij));
% iter = 0;
% lgd = {};
% h = 0;
% for icij = 1:length(cij) %[2,5] %1:length(cij)
%     iter = iter + 1;
%     for idep = 1:length(N_norm)
%         dep = depthlayer_C(idep);
%         [cij_model, dep] = ACFLNGBHE_to_cij(A_norm,C_norm,F_norm,L_norm,N_norm,Gc,Bc,Hc,Ec,depthlayer_C,dep);
%         cij_model = cij_model / 1e9;
% 
%         cij_resid = (cij_model-cij(icij).c) ./ cij_model;
%         cij_resid(isnan(cij_resid)) = 0;
%         cij_resid(2,1) = 0;
%         cij_resid(3,1) = 0;
%         cij_resid(3,2) = 0;
% 
%         misfit(idep) = sqrt(sum(cij_resid(:).^2))/9;
%     end
%     h(iter) = plot(misfit,-depthlayer_C,'-','color',clr(icij,:),'linewidth',4); hold on;
%     lgd{iter} = cij(icij).ref;
% end
% xlim([0 0.08]);
% ylabel('Depth');
% xlabel('C_{ij} Misfit');
% set(gca,'fontsize',16,'linewidth',2);
% legend(h,lgd,'location','southeastoutside','fontsize',12);
% 
% 
% % Radial Anisotropy
% xi = N_norm(1:length(N_norm))./L_norm(1:length(N_norm)); % N/L
%  
% figure(4); clf;
% plot(xi,-depthlayer_C,'k','linewidth',4); hold on;
% clr = lines(length(cij));
% iter = 0;
% lgd = {};
% h = 0;
% for icij = 1:length(cij) %[2,5] %1:length(cij)
%     iter = iter + 1;
%     xi_cij = cij_calc(icij).n/cij_calc(icij).l;
%     h(iter) = plot(xi_cij*ones(size(depthlayer_C)),-depthlayer_C,'-','color',clr(icij,:),'linewidth',4);
%     xlim([0.98 1.15]);
%     ylim([-20 -12]);
%     ylabel('Depth');
%     xlabel('\xi');
%     lgd{iter} = cij(icij).ref;
% end
% set(gca,'fontsize',16,'linewidth',2);
% legend(h,lgd,'location','southeast','fontsize',12);
% 
% save2pdf([figpath,'azi_depth_GBHE_bootstrap_alphF',num2str(alphF_GBH,'%2.1e'),'_nperm',num2str(nperm_S),'_nbs',num2str(nbs),'_ORALS_long_DAMP0_NORM_perc_MN86_GL_SPLIT_errprop_patch2','_cij_radial_compare.pdf'],4,1000)
% save2pdf([figpath,'azi_depth_GBHE_bootstrap_alphF',num2str(alphF_GBH,'%2.1e'),'_nperm',num2str(nperm_S),'_nbs',num2str(nbs),'_ORALS_long_DAMP0_NORM_perc_MN86_GL_SPLIT_errprop_patch2','_cij_matrix_compare.pdf'],5,1000)
% save2pdf([figpath,'azi_depth_GBHE_bootstrap_alphF',num2str(alphF_GBH,'%2.1e'),'_nperm',num2str(nperm_S),'_nbs',num2str(nbs),'_ORALS_long_DAMP0_NORM_perc_MN86_GL_SPLIT_errprop_patch2','_cij_misfit_compare.pdf'],6,1000)


%% plot inversion result for P2P strength and fastdir as function of Depth (SPLIT)
figure(103);clf;
clrs = lines(5);
fac = 0.7;

% STRENGTH
Nx = 2; Ny = 1;
sidegap = 0.10; topgap = 0.10; botgap = 0.10; vgap = 0.05; hgap = 0.05; cbar_bot = 0.04;
width = (1 - vgap*(Nx-1)-2*sidegap)/Nx; height = (1 - topgap - botgap - (Ny-1)*hgap)/Ny;
% set(gcf,'position',[ 118   290   1000   750]);
set(gcf,'position',[210     1   830   704]);
set(gcf,'color','w');
LBLFNT=18
ix = 1; iy = 1;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height*.75]);
hold on; set(gca,'fontsize',18,'linewidth',2);
% plot(strength_G_bs*100,depthlayer_G,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
h = patch([strength_B_l95*100 fliplr(strength_B_u95*100)],[depthlayer_B' fliplr(depthlayer_B')],clrs(2,:),'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([strength_B_l68*100 fliplr(strength_B_u68*100)],[depthlayer_B' fliplr(depthlayer_B')],clrs(2,:),'linestyle','none');
set(h, 'FaceAlpha', 0.45)
plot(strength_B_med*100,depthlayer_B,'-','color',clrs(2,:),'linewidth',4);hold on;
h = patch([strength_G_l95*100 fliplr(strength_G_u95*100)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([strength_G_l68*100 fliplr(strength_G_u68*100)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
h = patch([strength_H_l95*100 fliplr(strength_H_u95*100)],[depthlayer_H' fliplr(depthlayer_H')],clrs(3,:),'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([strength_H_l68*100 fliplr(strength_H_u68*100)],[depthlayer_H' fliplr(depthlayer_H')],clrs(3,:),'linestyle','none');
set(h, 'FaceAlpha', 0.45)
% plot(strength_C_bs*100,depthlayer_C,'-','color',[0.7 0.7 0.85],'linewidth',2); hold on;
h = patch([strength_C_l95*100 fliplr(strength_C_u95*100)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([strength_C_l68*100 fliplr(strength_C_u68*100)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
plot(strength_H_med*100,depthlayer_H,'-','color',clrs(3,:),'linewidth',4);hold on;
plot(strength_G_med*100,depthlayer_G,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);hold on;
plot(strength_C_med*100,depthlayer_C,'-','color',[0 0.4470 0.7410],'linewidth',4);hold on;
% plot(strength_G_bs(:,1)*100,depthlayer_G,'-k','linewidth',2);
xlim([0 8]); 
ylim([40 300]);
ylabel('Depth (km)','fontsize',18)
xlabel('Strength (%)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'YDir','reverse')

ix = 1.72; iy = 1;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
axes('position',[left,bot,width,height*.3]);
hold on; set(gca,'fontsize',18,'linewidth',2);
% plot(strength_G_bs*100,depthlayer_G,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
h = patch([strength_G_l95*100 fliplr(strength_G_u95*100)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([strength_G_l68*100 fliplr(strength_G_u68*100)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
h = patch([strength_B_l95*100 fliplr(strength_B_u95*100)],[depthlayer_B' fliplr(depthlayer_B')],clrs(2,:),'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([strength_B_l68*100 fliplr(strength_B_u68*100)],[depthlayer_B' fliplr(depthlayer_B')],clrs(2,:),'linestyle','none');
set(h, 'FaceAlpha', 0.45)
h = patch([strength_H_l95*100 fliplr(strength_H_u95*100)],[depthlayer_H' fliplr(depthlayer_H')],clrs(3,:),'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([strength_H_l68*100 fliplr(strength_H_u68*100)],[depthlayer_H' fliplr(depthlayer_H')],clrs(3,:),'linestyle','none');
set(h, 'FaceAlpha', 0.45)
% plot(strength_C_bs*100,depthlayer_C,'-','color',[0.7 0.7 0.85],'linewidth',2); hold on;
h = patch([strength_C_l95*100 fliplr(strength_C_u95*100)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([strength_C_l68*100 fliplr(strength_C_u68*100)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
plot(strength_B_med*100,depthlayer_B,'-','color',clrs(2,:),'linewidth',4);hold on;
plot(B_A_moho_7km*100,depth_moho_7km,'-.','color',clrs(2,:)*fac,'linewidth',4);hold on;
plot(strength_H_med*100,depthlayer_H,'-','color',clrs(3,:),'linewidth',4);hold on;
plot(strength_G_med*100,depthlayer_G,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);hold on;
plot(strength_C_med*100,depthlayer_C,'-','color',[0 0.4470 0.7410],'linewidth',4);hold on;
plot(E_N_moho_7km*100,depth_moho_7km,'-.','color',[0 0.4470 0.7410]*fac,'linewidth',4);hold on;
% plot(strength_G_bs(:,1)*100,depthlayer_G,'-k','linewidth',2);
clr = lines(length(cij_ind));
iter = 0;
bar_len = [15, 19]; %[15, 17, 19, 19];
for ic = cij_ind
    iter = iter+1;
    strength_g = cij_calc(ic).strength_g;
    strength_c = cij_calc(ic).strength_c;
%     h(ic) = plot([strength_g strength_g],[0 20],'-','color',clr(iter,:),'linewidth',4); hold on;
%     plot([strength_c strength_c],[0 20],'--','color',clr(iter,:),'linewidth',4);
    h(ic) = plot([strength_g strength_g],[0 bar_len(iter)],'-','color',[0.2660    0.6740    0.1880],'linewidth',4); hold on;
    plot([strength_c strength_c],[0 bar_len(iter)],'-','color',[0 0.4470 0.7410],'linewidth',4);
end
xlim([0 8]); 
ylim([11 40]);
% ylabel('Depth (km)','fontsize',18)
% xlabel('Strength (%)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'YDir','reverse','xticklabel',[])


% AZIMUTH
ix = 1;iy = 2;
left = sidegap + (iy-1)*(vgap+width);bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height*0.75]);
hold on;set(gca,'fontsize',18,'linewidth',2);hold on;set(gca,'fontsize',18,'linewidth',2);
% plot(fastdir_G_bs,depthlayer_G,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
h = patch([fastdir_G_l95 fliplr(fastdir_G_u95)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([fastdir_G_l68 fliplr(fastdir_G_u68)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
% h = patch([fastdir_B_l95 fliplr(fastdir_B_u95)],[depthlayer_B' fliplr(depthlayer_B')],clrs(2,:),'linestyle','none'); hold on;
% set(h, 'FaceAlpha', 0.2)
% h = patch([fastdir_B_l68 fliplr(fastdir_B_u68)],[depthlayer_B' fliplr(depthlayer_B')],clrs(2,:),'linestyle','none');
% set(h, 'FaceAlpha', 0.45)
% h = patch([fastdir_H_l95 fliplr(fastdir_H_u95)],[depthlayer_H' fliplr(depthlayer_H')],clrs(3,:),'linestyle','none'); hold on;
% set(h, 'FaceAlpha', 0.2)
% h = patch([fastdir_H_l68 fliplr(fastdir_H_u68)],[depthlayer_H' fliplr(depthlayer_H')],clrs(3,:),'linestyle','none');
% set(h, 'FaceAlpha', 0.45)
% plot(strength_C_bs*100,depthlayer_C,'-','color',[0.7 0.7 0.85],'linewidth',2); hold on;
h = patch([fastdir_C_l95 fliplr(fastdir_C_u95)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([fastdir_C_l68 fliplr(fastdir_C_u68)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
h1(2) = plot(fastdir_B_med, depthlayer_B,'-','color',clrs(2,:),'linewidth',4);hold on;
h1(3) = plot(fastdir_H_med, depthlayer_H,'-','color',clrs(3,:),'linewidth',4);hold on;
h1(1) = plot(fastdir_G_med, depthlayer_G,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);hold on;
h1(2) = plot(fastdir_C_med, depthlayer_C,'-','color',[0 0.4470 0.7410],'linewidth',4);hold on;
% plot(fastdir_G_bs(:,1),depthlayer_G,'-k','linewidth',2);
plot([ 78 78],[0 250],'k--','Linewidth',3);
APM = 296.74;
plot([ APM-180 APM-180 ],[0 400],'--','color',[0.6 0.6 0.6], 'Linewidth',3);
set(gca,'YDir','reverse')
xlim([45 150]);
ylim([40 300]);
xlabel('Azimuth (\circ)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
set(gca,'xtick',[0:30:180]);
% legend(h1,{'G/L (2\theta)','B/A (2\theta)','H/F (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');
% legend(h1,{'G/L (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');

ix = 1.72;iy = 2;
left = sidegap + (iy-1)*(vgap+width);bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height*0.3]);
hold on;set(gca,'fontsize',18,'linewidth',2);hold on;set(gca,'fontsize',18,'linewidth',2);
% plot(fastdir_G_bs,depthlayer_G,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
h = patch([fastdir_G_l95 fliplr(fastdir_G_u95)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([fastdir_G_l68 fliplr(fastdir_G_u68)],[depthlayer_G' fliplr(depthlayer_G')],[0.2660    0.6740    0.1880],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
% h = patch([fastdir_B_l95 fliplr(fastdir_B_u95)],[depthlayer_B' fliplr(depthlayer_B')],clrs(2,:),'linestyle','none'); hold on;
% set(h, 'FaceAlpha', 0.2)
% h = patch([fastdir_B_l68 fliplr(fastdir_B_u68)],[depthlayer_B' fliplr(depthlayer_B')],clrs(2,:),'linestyle','none');
% set(h, 'FaceAlpha', 0.45)
% h = patch([fastdir_H_l95 fliplr(fastdir_H_u95)],[depthlayer_H' fliplr(depthlayer_H')],clrs(3,:),'linestyle','none'); hold on;
% set(h, 'FaceAlpha', 0.2)
% h = patch([fastdir_H_l68 fliplr(fastdir_H_u68)],[depthlayer_H' fliplr(depthlayer_H')],clrs(3,:),'linestyle','none');
% set(h, 'FaceAlpha', 0.45)
% plot(strength_C_bs*100,depthlayer_C,'-','color',[0.7 0.7 0.85],'linewidth',2); hold on;
h = patch([fastdir_C_l95 fliplr(fastdir_C_u95)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none'); hold on;
set(h, 'FaceAlpha', 0.2)
h = patch([fastdir_C_l68 fliplr(fastdir_C_u68)],[depthlayer_C' fliplr(depthlayer_C')],[0 0.4470 0.7410],'linestyle','none');
set(h, 'FaceAlpha', 0.45)
h1(2) = plot(fastdir_B_med, depthlayer_B,'-','color',clrs(2,:),'linewidth',4);hold on;
plot(phi_B_moho_7km,depth_moho_7km,'-.','color',clrs(2,:)*fac,'linewidth',4);hold on;
h1(3) = plot(fastdir_H_med, depthlayer_H,'-','color',clrs(3,:),'linewidth',4);hold on;
h1(1) = plot(fastdir_G_med, depthlayer_G,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);hold on;
h1(2) = plot(fastdir_C_med, depthlayer_C,'-','color',[0 0.4470 0.7410],'linewidth',4);hold on;
plot(phi_E_moho_7km+90+45,depth_moho_7km,'-.','color',[0 0.4470 0.7410]*fac,'linewidth',4);hold on;
% plot(fastdir_G_bs(:,1),depthlayer_G,'-k','linewidth',2);
plot([ 78 78],[0 250],'k--','Linewidth',3);
plot([ 78 78]+45,[0 250],'k--','Linewidth',3);
APM = 296.74;
% plot([ APM-180 APM-180 ],[0 400],'--','color',[0.6 0.6 0.6], 'Linewidth',3);
set(gca,'YDir','reverse')
xlim([45 150]);
ylim([11 40]);
% xlabel('Azimuth (\circ)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
set(gca,'xticklabel',[]);
set(gca,'xtick',[0:30:180]);
% legend(h1,{'G/L (2\theta)','B/A (2\theta)','H/F (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');
legend(h1,{'G/L (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');



if isfig
    save2pdf([figpath,'BOOT_azi_depth_GBHE_H90_bootstrap_alphF',num2str(alphF_GBH,'%2.1e'),'_nbs',num2str(nbs),'_ORALS_long_DAMP0_NORM_perc_MN86_GL_SPLIT_errprop_patch.pdf'],103,1000)
    export_fig([figpath,'BOOT_azi_depth_GBHE_H90_bootstrap_alphF',num2str(alphF_GBH,'%2.1e'),'_nbs',num2str(nbs),'_ORALS_long_DAMP0_NORM_perc_MN86_GL_SPLIT_errprop_patch2.pdf'],'-pdf','-q100','-p0.02','-painters',103)
end


%% plot A and Phi fits
figure(104);clf;

% STRENGTH
Nx = 2; Ny = 1;
sidegap = 0.10; topgap = 0.10; botgap = 0.10; vgap = 0.05; hgap = 0.05; cbar_bot = 0.04;
width = (1 - vgap*(Nx-1)-2*sidegap)/Nx; height = (1 - topgap - botgap - (Ny-1)*hgap)/Ny;
%set(gcf,'position',[ 118   290   1000   750]);
set(gcf,'position',[210     1   830   704]);
set(gcf,'color','w');
LBLFNT=18
ix = 1; iy = 1;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height*.75]);
hold on; set(gca,'fontsize',18,'linewidth',2);
plot(new_A2_GBH_R0_bs*2*100,S0periods,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
plot(A2_R0_perc*2*100,S0periods,'o','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',12,'markerfacecolor',[0.2660    0.6740    0.1880]); hold on;
h = errorbar(A2_R0_perc*2*100,S0periods,err_A2_R0_perc*2*100/(err_pct),'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
plot(new_A2_GBH_R0_mean*2*100,S0periods,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);
ylim([min(S0periods)-10 max(S0periods)+10])
ylabel('Period (s)','fontsize',12);
xlabel('2A (%)','fontsize',12);
set(gca,'fontsize',12,'linewidth',1.5,'ydir','reverse')
xlim([0 5])
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'YDir','reverse')

ix = 1.72; iy = 1;
left = sidegap + (iy-1)*(vgap+width); bot = botgap + (ix-1)*(hgap+height);
axes('position',[left,bot,width,height*.3]);
hold on; set(gca,'fontsize',18,'linewidth',2);
plot(new_A2_GBH_R1_bs*2*100,S1periods-.1,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
plot(new_A2_G_L_bs*2*100,T0periods+.1,'-','color',[0.7 0.85 0.7],'linewidth',2)
plot(new_A4_C_L_bs*2*100,T0periods,'-','color',[0.6824 0.8353 0.9765],'linewidth',2)
plot(A2_R1_perc*2*100,S1periods-.1,'o','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',12,'markerfacecolor',[0.2660    0.6740    0.1880]); hold on;
plot(A2_L_perc*2*100,T0periods+.1,'^','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',12,'markerfacecolor',[0.2660    0.6740    0.1880]);
plot(A4_L_perc*2*100,T0periods,'^','color',[0 0.4470 0.7410],'linewidth',2,'markersize',12,'markerfacecolor',[0 0.4470 0.7410]);
h = errorbar(A2_R1_perc*2*100,S1periods-.1,err_A2_R1_perc*2*100/(err_pct),'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(A2_L_perc*2*100,T0periods+.1,err_A2_L_perc*2*100/(err_pct),'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(A4_L_perc*2*100,T0periods,err_A4_L_perc*2*100/(err_pct),'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
plot(new_A2_GBH_R1_mean*2*100,S1periods-.1,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);
plot(new_A2_G_L_mean*2*100,T0periods+.1,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);
plot(new_A4_C_L_mean*2*100,T0periods,'-','color',[0 0.4470 0.7410],'linewidth',4);
ylim([min(S1periods)-.25 max(S1periods)+.25])
% ylabel('Period (s)','fontsize',12);
% xlabel('2A (%)','fontsize',12);
set(gca,'xticklabel',[],'fontsize',12,'linewidth',1.5,'ydir','reverse')
xlim([0 5])
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'YDir','reverse','xticklabel',[])



% AZIMUTH
ix = 1;iy = 2;
left = sidegap + (iy-1)*(vgap+width);bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height*0.75]);
hold on;set(gca,'fontsize',18,'linewidth',2);hold on;set(gca,'fontsize',18,'linewidth',2);
plot(new_phi2_GBH_R0_bs,S0periods,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
plot(phi2_R0,S0periods,'o','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',12,'markerfacecolor',[0.2660    0.6740    0.1880]);hold on;
h = errorbar(phi2_R0,S0periods,err_phi2_R0,'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
plot(new_phi2_GBH_R0_mean,S0periods,'-','color',[0.2660    0.6740    0.1880],'linewidth',4)
xlim([50 200]);
ylim([min(S0periods)-10 max(S0periods)+10])
% ylabel('Period (s)','fontsize',12);
xlabel('\psi (\circ)','fontsize',12)
plot([ 78 78],[0 250],'k--','Linewidth',3);
APM = 296.74;
plot([ APM-180 APM-180 ],[0 400],'--','color',[0.6 0.6 0.6], 'Linewidth',3);
set(gca,'YDir','reverse')
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
set(gca,'xtick',[0:30:180]);
% legend(h1,{'G/L (2\theta)','B/A (2\theta)','H/F (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');
% legend(h1,{'G/L (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');

ix = 1.72;iy = 2;
left = sidegap + (iy-1)*(vgap+width);bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height*0.3]);
hold on;set(gca,'fontsize',18,'linewidth',2);hold on;set(gca,'fontsize',18,'linewidth',2);
plot(new_phi2_GBH_R1_bs,S1periods,'-','color',[0.7 0.85 0.7],'linewidth',2); hold on;
plot(new_phi2_G_L_bs,T0periods,'-','color',[0.7 0.85 0.7],'linewidth',2)
plot(new_phi4_C_L_bs,T0periods,'-','color',[0.6824 0.8353 0.9765],'linewidth',2)
plot(phi2_R1,S1periods,'o','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',12,'markerfacecolor',[0.2660    0.6740    0.1880]);hold on;
plot(phi2_L,T0periods,'^','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',12,'markerfacecolor',[0.2660    0.6740    0.1880]);hold on;
plot(phi4_L2,T0periods,'^','color',[0 0.4470 0.7410],'linewidth',2,'markersize',12,'markerfacecolor',[0 0.4470 0.7410]);hold on;
h = errorbar(phi2_R1,S1periods,err_phi2_R1,'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(phi2_L,T0periods,err_phi2_L,'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(phi4_L2,T0periods,err_phi4_L,'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
plot(new_phi2_GBH_R1_mean,S1periods,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);
plot(new_phi2_G_L_mean,T0periods,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);
plot(new_phi4_C_L_mean,T0periods,'-','color',[0 0.4470 0.7410],'linewidth',4);
xlim([50 200]);
ylim([min(S1periods)-0.25 max(S1periods)+0.25])
% ylabel('Period (s)','fontsize',12);
plot([ 78 78],[0 250],'k--','Linewidth',3);
plot([ 78 78]+45,[0 250],'k--','Linewidth',3);
plot([ 78 78]+90,[0 250],'k--','Linewidth',3);
APM = 296.74;
% plot([ APM-180 APM-180 ],[0 400],'--','color',[0.6 0.6 0.6], 'Linewidth',3);
set(gca,'YDir','reverse')
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
set(gca,'xticklabel',[]);
% legend(h1,{'G/L (2\theta)','B/A (2\theta)','H/F (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');
% legend(h1,{'G/L (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off');

if isfig
    save2pdf([figpath,'BOOT_azi_depth_GBHE_H90_bootstrap_alphF',num2str(alphF_GBH,'%2.1e'),'_nbs',num2str(nbs),'_ORALS_long_DAMP0_NORM_perc_MN86_A_PHI_errprop_patch.pdf'],104,1000)
    export_fig([figpath,'BOOT_azi_depth_GBHE_H90_bootstrap_alphF',num2str(alphF_GBH,'%2.1e'),'_nbs',num2str(nbs),'_ORALS_long_DAMP0_NORM_perc_MN86_A_PHI_errprop_patch2.pdf'],'-pdf','-q100','-p0.02','-painters',104)
end