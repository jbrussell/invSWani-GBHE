clear;
close all;
%% invert for
% (Gc, Gs); (Bc, Bs); (Hc, Hs); (Ec, Es)
%
% NOTE:
% In this code the elastic parameter "E" for 4-theta anisotropy is called "C" to 
% match Montagner & Nataf (1986). But it's confusing, and most literature 
% uses "E".
%
%% Parameters
% parameter_ACFLN
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
phi_patty = 78; % fossil spreading direction

%%
%%%%%%% INVERSION PARAMETERS %%%%%%%
model_depth_G = 400;
model_depth_B = 400;
model_depth_H = model_depth_B;
model_depth_C = 40; %25;

% smoothing (second derivative)
    alphF_GBH = 9e2; %9e2; %2e3 
    alphF_C = 3e3; %9e2; %3e3; 
    
% Flatness (first derivative);
    alphJ_C = 2e3; %2e3; %3e3; 
    
% Damp below certain depth
    G_DAMP = 5e1; %5e1; %3e2;
    dep_zero_damp = 300;

% Norm Damping
    alphH_GBH = 1e0; %0.005 %good : 0.001 %smaller fit more
    alphH_C = 5e1; %5e1; %0.005 %good : 0.001 %smaller fit more
    
% Scaling Ratios
    epsilonBG = 1e3; %1e3; %1e0;  % Enforce B/G ratio
    BGratio = 1.25; %1.25;
    epsilonHG = 1e3; %1e3; %1e0;  % Enforce H/G ratio
    HGratio = -0.11; %0.11; %0.25; %0.25; %1.25;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isfig = 0;
is_RMS = 1;
issavemat = 0;

%% get G matrix

T0_Gmat = load(Gmat_T0path);
S0_Gmat = load(Gmat_S0path);
S1_Gmat = load(Gmat_S1path);

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
% vpv_cutoff = 8000; % 4000 (crust)
vpv_cutoff = 7500; % 4000 (crust)
I = find(Velmodel.card.vpv > vpv_cutoff);
qeqeqe = length(Velmodel.card.vpv)-I(end);
if vpv_cutoff == 1500
    top = 'noh20';
elseif vpv_cutoff == 3000
    top = 'noh20sed';
elseif vpv_cutoff == 8000
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
GGG2_R0 = midpts(S0_Gmat.Gmatrix.L(:,IdepG)).*drG'.*(c_S0'./U_S0'.*rho'.*Vsv'.^2);
GGG2_R1 = midpts(S1_Gmat.Gmatrix.L(:,IdepG)).*drG'.*(c_S1'./U_S1'.*rho'.*Vsv'.^2);
GGG2_L = midpts(-T0_Gmat.Gmatrix.L(:,IdepG)).*drG'.*(c_T0'./U_T0'.*rho'.*Vsv'.^2);

% B (A) Rayleigh
BBB2_R0 = midpts(S0_Gmat.Gmatrix.A(:,IdepB)).*drB'.*(c_S0'./U_S0'.*rho'.*Vph'.^2);
BBB2_R1 = midpts(S1_Gmat.Gmatrix.A(:,IdepB)).*drB'.*(c_S1'./U_S1'.*rho'.*Vph'.^2);

% H (F) Rayleigh
HHH2_R0 = midpts(S0_Gmat.Gmatrix.F(:,IdepH)).*drH'.*(c_S0'./U_S0'.*rho'.*eta'.*(Vph'.^2 - 2*Vsv'.^2));
HHH2_R1 = midpts(S1_Gmat.Gmatrix.F(:,IdepH)).*drH'.*(c_S1'./U_S1'.*rho'.*eta'.*(Vph'.^2 - 2*Vsv'.^2));

% 4 theta
% C (A) Rayleigh, (-N) Love
CCC4 = midpts(-T0_Gmat.Gmatrix.N(:,IdepC)).*drC'.*(c_T0'./U_T0'.*rho_C'.*Vsh_C'.^2);

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
end
GGG2 = [GGG2_R0;
        GGG2_R1;
        GGG2_L];

%% Build data vectors
I_R0 = aniso_R.mode_br == 0;
I_R1 = aniso_R.mode_br == 1;

% 2 theta - Cosine part
c2_R0 = c2_R(I_R0);
c2_R1 = c2_R(I_R1);
    c2_R = [c2_R0, c2_R1]';
c2_L = c2_L';
c2_std_R0 = c2_std_R(I_R0);
c2_std_R1 = c2_std_R(I_R1);
    c2_std_R = [c2_std_R0, c2_std_R1]';
c2_std_L = c2_std_L';
% 2 theta - Sine part
s2_R0 = s2_R(I_R0);
s2_R1 = s2_R(I_R1);
    s2_R = [s2_R0, s2_R1]';
s2_L = s2_L';
s2_std_R0 = s2_std_R(I_R0);
s2_std_R1 = s2_std_R(I_R1);
    s2_std_R = [s2_std_R0, s2_std_R1]';
s2_std_L = s2_std_L';
% 4 theta - Cosine part
c4_R = c4_R';
c4_L = c4_L';
c4_std_R = c4_std_R';
c4_std_L = c4_std_L';
% 4 theta - Sine part
s4_R = s4_R';
s4_L = s4_L';
s4_std_R = s4_std_R';
s4_std_L = s4_std_L';

%% INVERSION

% Setup damping to zero beneath specified depth (JOSH)
GBH_EYE = eye(nlayers_G+nlayers_B+nlayers_H);
I_DAMP = find(midpts(S0_Gmat.depthlayer(IdepG)')' >= dep_zero_damp);
G_EYE_DAMP = [GBH_EYE(I_DAMP,:)];

nlayer_G_DAMP = size(G_EYE_DAMP,1);

% Damp linearly from 0 to G_DAMP
G_DAMP_vec = linspace(0,G_DAMP,nlayer_G_DAMP)';


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

% Cosine
d_GBH_c2 = [c2_R;
            c2_L;
           zeros(nlayer_G_DAMP,1);
           zeros(nlayers_B,1); % B/G ratio
           zeros(nlayers_H,1)]; % H/G ratio

GGGBBBHHH2 = [GGG2_R0, BBB2_R0, HHH2_R0;
              GGG2_R1, BBB2_R1, HHH2_R1;
              GGG2_L, zeros(size(GGG2_L)), zeros(size(GGG2_L));
            G_EYE_DAMP;
           -BGratio*eye(nlayers_G), eye(nlayers_B), zeros(nlayers_H); % B/G ratio
           -HGratio*eye(nlayers_G), zeros(nlayers_B), eye(nlayers_H)]; % H/G ratio

% invert for MGc
std_GBH_c2 = [c2_std_R;
              c2_std_L;
             ones(nlayer_G_DAMP,1)./G_DAMP_vec;
             ones(nlayers_G,1)/epsilonBG; % B/G ratio
             ones(nlayers_G,1)/epsilonHG]; % H/G ratio
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

new_c2_G_R0 = GGG2_R0*MG_c2;
new_c2_B_R0 = BBB2_R0*MB_c2;
new_c2_H_R0 = HHH2_R0*MH_c2;
new_c2_GBH_R0 = new_c2_G_R0 + new_c2_B_R0 + new_c2_H_R0;
new_c2_GBH_R0_test = GGGBBBHHH2*MGBH_c2;
res_c2_GBH_R0 = c2_R0' - new_c2_GBH_R0;

new_c2_G_R1 = GGG2_R1*MG_c2;
new_c2_B_R1 = BBB2_R1*MB_c2;
new_c2_H_R1 = HHH2_R1*MH_c2;
new_c2_GBH_R1 = new_c2_G_R1 + new_c2_B_R1 + new_c2_H_R1;
new_c2_GBH_R1_test = GGGBBBHHH2*MGBH_c2;
res_c2_GBH_R1 = c2_R1' - new_c2_GBH_R1;

new_c2_G_L = GGG2_L*MG_c2;
res_c2_G_L = c2_L' - new_c2_G_L;




% Sine
d_GBH_s2 = [s2_R;
            s2_L;
           zeros(nlayer_G_DAMP,1);
           zeros(nlayers_B,1); % B/G ratio
           zeros(nlayers_H,1)]; % H/G ratio

% invert for MGs
std_GBH_s2 = [s2_std_R;
              s2_std_L;
             ones(nlayer_G_DAMP,1)./G_DAMP_vec;
             ones(nlayers_G,1)/epsilonBG; % B/G ratio
             ones(nlayers_G,1)/epsilonHG]; % H/G ratio
C_std_GBH_s2 = diag(std_GBH_s2.^2) ; % std matrix
GGBBHH2 = inv(GGGBBBHHH2'* inv(C_std_GBH_s2)*GGGBBBHHH2+ Q_GBH) * GGGBBBHHH2';
MGBH_s2 = GGBBHH2*inv(C_std_GBH_s2)*d_GBH_s2;

MG_s2 = MGBH_s2(1:nlayers_G);
MB_s2 = MGBH_s2(nlayers_G+1:nlayers_G+nlayers_B);
MH_s2 = MGBH_s2(nlayers_G+nlayers_B+1:nlayers_G+nlayers_B+nlayers_H);
% display(MB_s2./MG_s2);
% display(MH_s2./MG_s2);

new_s2_G_R0 = GGG2_R0*MG_s2;
new_s2_B_R0 = BBB2_R0*MB_s2;
new_s2_H_R0 = HHH2_R0*MH_s2;
new_s2_GBH_R0 = new_s2_G_R0 + new_s2_B_R0 + new_s2_H_R0;
new_s2_GBH_R0_test = GGGBBBHHH2*MGBH_s2;
res_s2_GBH_R0 = s2_R0' - new_s2_GBH_R0;

new_s2_G_R1 = GGG2_R1*MG_s2;
new_s2_B_R1 = BBB2_R1*MB_s2;
new_s2_H_R1 = HHH2_R1*MH_s2;
new_s2_GBH_R1 = new_s2_G_R1 + new_s2_B_R1 + new_s2_H_R1;
new_s2_GBH_R1_test = GGGBBBHHH2*MGBH_s2;
res_s2_GBH_R1 = s2_R1' - new_s2_GBH_R1;

new_s2_G_L = GGG2_L*MG_s2;
res_s2_G_L = s2_L_save' - new_s2_G_L;

%% C

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
d_C_c4 = [c4_L];

% invert for MCc
std_C_c4 = [c4_std_L];
C_std_C_c4 = diag(std_C_c4.^2) ; % std matrix
H0_C = H0./norm(H0).*norm(CCC4');
F0_C = F0./norm(F0).*norm(CCC4');
J0_C = J0./norm(J0).*norm(CCC4');
H_C = alphH_C * H0_C;
F_C = alphF_C * F0_C;
J_C = alphJ_C * J0_C;
Q_C = F_C'*F_C + H_C'*H_C + J_C'*J_C;
CC4 = inv(CCC4'* inv(C_std_C_c4)*CCC4+ Q_C) * CCC4';
MC_c4 = CC4*inv(C_std_C_c4)*d_C_c4;

new_c4_C_L = CCC4*MC_c4;
res_c4_C_L = c4_L' - new_c4_C_L;

% Sine
d_C_s4 = [s4_L];

% invert for MCs
std_C_s4 = [s4_std_L];
C_std_C_s4 = diag(std_C_s4.^2) ; % std matrix
CC4 = inv(CCC4'* inv(C_std_C_s4)*CCC4 + Q_C) * CCC4';
MC_s4 = CC4*inv(C_std_C_s4)*d_C_s4;

new_s4_C_L = CCC4*MC_s4;
res_s4_C_L = s4_L_save' - new_s4_C_L;

%% get strength and direction from MGs and MGc =================

% G
strength_G = sqrt(MG_c2.^2+MG_s2.^2); % 2*dV/V = G/L
%strength = 2*strength ./ normalize_strength(qeqeqe:end);  
% strength_G = strength_G ;%./ L_norm;
fastdir_G = 0.5*atan2d(MG_s2,MG_c2);
fastdir_vec = [];
for idep = 1:length(depthlayer_G)
    
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
    
    fastdir_vec(1) = fastdir_H(idep);
    fastdir_vec(2) = fastdir_H(idep)+180;
    fastdir_vec(3) = fastdir_H(idep)-180;
    [~, I] = min(abs(fastdir_vec-phi_patty-90));
    fastdir2_H(idep) = fastdir_vec(I);
       
end
ind = find(fastdir2_H<0);
fastdir2_H(ind) = fastdir2_H(ind)+180;

% C
strength_C = sqrt(MC_c4.^2+MC_s4.^2);
fastdir_C = (1/4)*atan2d(-MC_s4,-MC_c4);
ind = find(fastdir_C<0);
fastdir_C(ind) = fastdir_C(ind)+180;
for idep = 1:length(depthlayer_C)
    
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

%% plot inversion result for P2P strength and fastdir as function of Depth (SPLIT)
figure(103);clf;

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
h1(2) = plot(strength_B*100,depthlayer_B,'-','color',[1 0.7 0.7],'linewidth',4);hold on;
h1(3) = plot(strength_H*100,depthlayer_H,'-','color',[0.7 0.7 1],'linewidth',4);hold on;
h1(1) = plot(strength_G*100,depthlayer_G,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);hold on;
h1(4) = plot(strength_C*100,depthlayer_C,'-','color',[0 0.4470 0.7410],'linewidth',4);hold on;
% plot(strength_G_bs(:,1)*100,depthlayer_G,'-k','linewidth',2);
xlim([0 8]); 
ylim([0 300]);
ylabel('Depth (km)','fontsize',18)
xlabel('Strength (%)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
% title('$G/L = 2 \,\delta V_{SV}/V_{SV}$','interpreter','Latex','fontsize',30)
set(gca,'YDir','reverse')
legend(h1,{'G/L (2\theta)','B/A (2\theta)','H/F (2\theta)','E/N (4\theta)'},'fontsize',18,'box','off','Location','southeast');



% AZIMUTH
ix = 1;iy = 2;
left = sidegap + (iy-1)*(vgap+width);bot = botgap + (ix-1)*(hgap+height);
subplot('position',[left,bot,width,height]);
hold on;set(gca,'fontsize',18,'linewidth',2);hold on;set(gca,'fontsize',18,'linewidth',2);
h1(2) = plot(fastdir2_B, depthlayer_B,'-','color',[1 0.7 0.7],'linewidth',4);hold on;
h1(3) = plot(fastdir2_H, depthlayer_H,'-','color',[0.7 0.7 1],'linewidth',4);hold on;
h1(1) = plot(fastdir2_G, depthlayer_G,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);hold on;
h1(4) = plot(fastdir2_C, depthlayer_C,'-','color',[0 0.4470 0.7410],'linewidth',4);hold on;
% plot(fastdir_G_bs(:,1),depthlayer_G,'-k','linewidth',2);
plot([ 78 78],[0 250],'k--','Linewidth',3);
APM = 296.74;
plot([ APM-180 APM-180 ],[0 400],'--','color',[0.6 0.6 0.6], 'Linewidth',3);
set(gca,'YDir','reverse')
xlim([45 150]);
ylim([0 300]);
xlabel('Azimuth (\circ)','fontsize',18)
set(gca,'box','on','LineWidth', 2,'xminortick','on','yminortick','on','fontsize', LBLFNT);
set(gca,'xtick',[0:30:180]);

if isfig
    save2pdf([figpath,'azi_depth_GBHE.pdf'],103,1000)
%     export_fig([figpath,'azi_depth_GBHE_2.pdf'],'-pdf','-q100','-p0.02','-painters',103)
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
plot(A2_R0_perc*2*100,S0periods,'o','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',12,'markerfacecolor',[0.2660    0.6740    0.1880]); hold on;
h = errorbar(A2_R0_perc*2*100,S0periods,err_A2_R0_perc*2*100/(err_pct),'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
plot(new_A2_GBH_R0*2*100,S0periods,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);
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
plot(A2_R1_perc*2*100,S1periods-.1,'o','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',12,'markerfacecolor',[0.2660    0.6740    0.1880]); hold on;
plot(A2_L_perc*2*100,T0periods+.1,'^','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',12,'markerfacecolor',[0.2660    0.6740    0.1880]);
plot(A4_L_perc*2*100,T0periods,'^','color',[0 0.4470 0.7410],'linewidth',2,'markersize',12,'markerfacecolor',[0 0.4470 0.7410]);
h = errorbar(A2_R1_perc*2*100,S1periods-.1,err_A2_R1_perc*2*100/(err_pct),'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(A2_L_perc*2*100,T0periods+.1,err_A2_L_perc*2*100/(err_pct),'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(A4_L_perc*2*100,T0periods,err_A4_L_perc*2*100/(err_pct),'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
plot(new_A2_GBH_R1*2*100,S1periods-.1,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);
plot(new_A2_G_L*2*100,T0periods+.1,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);
plot(new_A4_C_L*2*100,T0periods,'-','color',[0 0.4470 0.7410],'linewidth',4);
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
plot(phi2_R0,S0periods,'o','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',12,'markerfacecolor',[0.2660    0.6740    0.1880]);hold on;
h = errorbar(phi2_R0,S0periods,err_phi2_R0,'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
plot(new_phi2_GBH_R0,S0periods,'-','color',[0.2660    0.6740    0.1880],'linewidth',4)
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
plot(phi2_R1,S1periods,'o','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',12,'markerfacecolor',[0.2660    0.6740    0.1880]);hold on;
plot(phi2_L,T0periods,'^','color',[0.2660    0.6740    0.1880],'linewidth',2,'markersize',12,'markerfacecolor',[0.2660    0.6740    0.1880]);hold on;
plot(phi4_L2,T0periods,'^','color',[0 0.4470 0.7410],'linewidth',2,'markersize',12,'markerfacecolor',[0 0.4470 0.7410]);hold on;
h = errorbar(phi2_R1,S1periods,err_phi2_R1,'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(phi2_L,T0periods,err_phi2_L,'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
h = errorbar(phi4_L2,T0periods,err_phi4_L,'horizontal','.','linewidth',2,'color','k');
h.CapSize = 0;
plot(new_phi2_GBH_R1,S1periods,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);
plot(new_phi2_G_L,T0periods,'-','color',[0.2660    0.6740    0.1880],'linewidth',4);
plot(new_phi4_C_L2,T0periods,'-','color',[0 0.4470 0.7410],'linewidth',4);
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
    save2pdf([figpath,'BOOT_azi_depth_GBHE_bootstrap_alphF',num2str(alphF_GBH,'%2.1e'),'_nbs',num2str(nbs),'_ORALS_long_DAMP0_NORM_perc_MN86_A_PHI_errprop_patch.pdf'],104,1000)
%     export_fig([figpath,'BOOT_azi_depth_GBHE_bootstrap_alphF',num2str(alphF_GBH,'%2.1e'),'_nbs',num2str(nbs),'_ORALS_long_DAMP0_NORM_perc_MN86_A_PHI_errprop_patch2.pdf'],'-pdf','-q100','-p0.02','-painters',104)
end

%% Save matfile of model
if issavemat
    matpath = ['./mats_singlestep/'];
    if ~exist(matpath)
        mkdir(matpath);
    end

    model.depthlayer_G = depthlayer_G;
    model.strength_G = strength_G_med;
    model.fastdir_G = fastdir_G_med;
    model.depthlayer_B = depthlayer_B;
    model.strength_B = strength_B_med;
    model.fastdir_B = fastdir_B_med;
    model.depthlayer_H = depthlayer_H;
    model.strength_H = strength_H_med;
    model.fastdir_H = fastdir_H_med;
    model.depthlayer_E = depthlayer_C;
    model.strength_E = strength_C_med;
    model.fastdir_E = fastdir_C_med;

    model.forward.S0periods = S0periods;
    model.forward.A2_GBH_R0 = new_A2_GBH_R0_mean;
    model.forward.phi2_GBH_R0 = new_phi2_GBH_R0_mean;
    model.forward.S1periods = S1periods;
    model.forward.A2_GBH_R1 = new_A2_GBH_R1_mean;
    model.forward.phi2_GBH_R1 = new_phi2_GBH_R1_mean;
    model.forward.T0periods = T0periods;
    model.forward.A2_G_L = new_A2_G_L_mean;
    model.forward.phi2_G_L = new_phi2_G_L_mean;
    model.forward.A4_E_L = new_A4_C_L_mean;
    model.forward.phi4_E_L = new_phi4_C_L_mean;

    model.data.S0periods = S0periods;
    model.data.A2_R0_perc = A2_R0_perc;
    model.data.err_A2_R0_perc = err_A2_R0_perc;
    model.data.phi2_R0 = phi2_R0;
    model.data.err_phi2_R0 = err_phi2_R0;
    model.data.S1periods = S1periods;
    model.data.A2_R1_perc = A2_R1_perc;
    model.data.err_A2_R1_perc = err_A2_R1_perc;
    model.data.phi2_R1 = phi2_R1;
    model.data.err_phi2_R1 = err_phi2_R1;
    model.data.T0periods = T0periods;
    model.data.A2_L_perc = A2_L_perc;
    model.data.err_A2_L_perc = err_A2_L_perc;
    model.data.phi2_L = phi2_L;
    model.data.err_phi2_L = err_phi2_L;
    model.data.A4_L_perc = A4_L_perc;
    model.data.err_A4_L_perc = err_A4_L_perc;
    model.data.phi4_L = phi4_L2;
    model.data.err_phi4_L = err_phi4_L;

    save([matpath,param.CARDID,'.mat'],'model');
end
