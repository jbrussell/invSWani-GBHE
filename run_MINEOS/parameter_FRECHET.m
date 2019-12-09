%% Setup Parameters for running MINEOS to calculate senstivity kernels, dispersion, and synthetics

%clear all;
% prepend functions directory to MATLAB path
fullMAINpath = mfilename('fullpath');
functionspath = [fullMAINpath(1:regexp(fullMAINpath,mfilename)-1),'functions'];
addpath(functionspath);

path2runMINEOS = './'; % Path to this folder
path2BIN = '../FORTRAN/bin'; % Path to fortran binaries

% Data for azimuthal anisotropy inversion
% param.DATAmat_path = ['./DATA/azi_measurements_RL_5_150s.mat'];
% param.DATAmat_path = ['./DATA/azi_measurements_RL_5_150s_test.mat']; % jbr 12/19
% param.DATAmat_path = ['./DATA/data_AmbGSDFLoveRayl_AGU19_5_150s.mat'];
% param.DATAmat_path = ['./DATA/data_AmbGSDFLoveRayl_AGU19_5_150s_loveT1.mat'];
% param.DATAmat_path = ['./DATA/data_AmbGSDFLoveRayl_AGU19_5_150s_oldAmbS0.mat'];
param.DATAmat_path = ['./DATA/data_AmbGSDFLoveRayl_AGU19_5_150s_oldAmbS0_loveT1.mat'];

% Mineos table parameters
maxN = 400000; % Estimate of max number of modes
minF = 0;
maxF = 200.05; % max frequency in mHz; %10.1; %250.05; %333.4; %500.05; %200.05; %%150.05; %50.05;
minL = 0;
maxL = 50000;
N_modes = 2; % <0 uses all mode branches, 1=fundamental only -------- JOSH 8/22/15
% param.CARDID = 'Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19'; %'Nomelt_taper_eta_crust_INV_noQ'; %'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays_layer2_5_150s_goodkerns3';
param.CARDID = 'Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19_ORCAiso_INV';

% (1 => yes, 0 => no)
SONLY = 0; %Spheroidal modes? (RAYLEIGH)
TONLY = 1; %Toroidal modes? (LOVE)
branch = 0;

% % for plotting kernels
% % param.periods = round(logspace(log10(5),log10(200),15));

% Do this in the script instead...
% (1 => yes, 0 => no)
% SONLY = 0; %Spheroidal modes? (RAYLEIGH)
% TONLY = 1; %Toroidal modes? (LOVE)
% branch = 0;
if exist('mode.txt') == 2
    fid = fopen('mode.txt');
    temp = textscan(fid,'%.0f %.0f %.0f');
    fclose(fid);
    SONLY = temp{1};
    TONLY = temp{2};
    branch = temp{3};
end

ch_mode = 0; % (DO NOT CHANGE) mode branch to check for missed eigenfrequencies 0 => T0 ------- JOSH 10/7/15

%% Parameters for idagrn synthetics
LENGTH_HR = 1.0; %1.0; % length of seismogram in hours
DT = 1.0; % 1/samplerate
eventfile = 'evt_201404131236';
stationfile = 'stations.stn';


%%
% Setup idagrn paths
param.IDAGRN = [path2runMINEOS,'/IDAGRN/'];
param.EVTPATH = [param.IDAGRN,'EVT_FILES/',eventfile];
param.STAPATH = [param.IDAGRN,'STATION/',stationfile];
param.SYNTH_OUT = [param.IDAGRN,'SYNTH/',param.CARDID,'_b',num2str(N_modes),'/',eventfile,'/'];
if ~exist(param.SYNTH_OUT)
    mkdir(param.SYNTH_OUT);
end

%%
if SONLY == 1 && TONLY == 0
    param.TYPE = 'S';
elseif SONLY == 0 && TONLY == 1
    param.TYPE = 'T';
else
    error('Choose SONLY or TONLY, not both');
    
end

% Setup Parameters for Initial Model
param.CARD = [param.CARDID,'.card'];
param.CARDPATH  = [path2runMINEOS,'/CARDS/'];
param.TABLEPATH = [path2runMINEOS,'/MODE/TABLES/'];
param.MODEPATH  = [path2runMINEOS,'/MODE/TABLES/MODE.in/'];
if ~exist(param.MODEPATH)
    mkdir(param.MODEPATH);
end
param.RUNPATH = pwd;

%% create dir for output MINEOS automatically, doesn't need to be changed.
CARDTABLE = [param.TABLEPATH,param.CARDID,'/tables/'];
if ~exist(CARDTABLE)
    mkdir([param.TABLEPATH,param.CARDID])
    mkdir(CARDTABLE)
end

%% setup Parameters for kernels
param.frechetpath = [path2runMINEOS,'/MODE/FRECHET/',param.CARDID,'/'];
param.frechet = [path2runMINEOS,'/MODE/FRECHET/'];
if ~exist(param.frechetpath) 
    mkdir(param.frechetpath)
end

%% setup Parameters for eigenfunctions
param.eigpath = [path2runMINEOS,'/MODE/EIGEN/',param.CARDID,'/'];

if ~exist(param.eigpath) 
    mkdir(param.eigpath)
end

%% setup Parameters for Dispersion
param.disperspath = [path2runMINEOS,'/MODE/DISPERSION/',param.CARDID,'/'];

if ~exist(param.disperspath) 
    mkdir(param.disperspath)
end

%% Turn on if only want to calculate S or T or both for mineous
param.SMODEIN = ['s.mode',num2str(floor(minF)),'_',num2str(floor(maxF)),'_b',num2str(N_modes)];
param.STYPEID = ['s',num2str(floor(minF)),'to',num2str(floor(maxF))];
param.TMODEIN = ['t.mode',num2str(floor(minF)),'_',num2str(floor(maxF)),'_b',num2str(N_modes)];
param.TTYPEID = ['t',num2str(floor(minF)),'to',num2str(floor(maxF))];%'t0to150';

%% Setup paths to FORTRAN binaries
PATH = getenv('PATH');
if isempty(strfind(PATH,path2BIN))
%     setenv('PATH', [PATH,':',path2BIN]);
    setenv('PATH', [path2BIN,':',PATH]);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%  START SETUP FOR AZIMUTHAL ANISOTROPY %%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Make Card Model
param.modeldirpath = ['./CARD_MODELS/'];
if ~exist(param.modeldirpath) 
    mkdir(param.modeldirpath)
end

param.modelpath = [param.modeldirpath,param.CARDID,'.mat'];
if ~exist(param.modelpath)
    card = read_model_card([param.CARDPATH,param.CARD]);
    save(param.modelpath,'card');
end

%% Load Anisotropy Measurements
% snr_tol = 0; % 5
% r_tol = 200; %200; % 100 km
% err_tol = 100; %0.7; %100;

% Load periods
load(param.DATAmat_path);
param.S0periods = data.rayl.periods_ani(data.rayl.mode_br_ani == 0);
param.S1periods = data.rayl.periods_ani(data.rayl.mode_br_ani == 1);
param.T0periods = data.love.periods_ani(data.love.mode_br_ani == 0);
param.T1periods = data.love.periods_ani(data.love.mode_br_ani == 1);

%% Load G matrix
% Sbranch = 1; % Branch Number
% Tbranch = 0; % Branch Number
% Spheroidal
param.Gmat_S0path = [param.frechetpath,'Gmatrix_S0_',num2str(param.S0periods(1)),'_',num2str(param.S0periods(end)),'s_',param.CARDID,'.mat'];
param.Gmat_S1path = [param.frechetpath,'Gmatrix_S1_',num2str(param.S1periods(1)),'_',num2str(param.S1periods(end)),'s_',param.CARDID,'.mat'];

% Toroidal
param.Gmat_T0path = [param.frechetpath,'Gmatrix_T0_',num2str(param.T0periods(1)),'_',num2str(param.T0periods(end)),'s_',param.CARDID,'.mat'];
if ~isempty(param.T1periods)
    param.Gmat_T1path = [param.frechetpath,'Gmatrix_T1_',num2str(param.T1periods(1)),'_',num2str(param.T1periods(end)),'s_',param.CARDID,'.mat'];
end
%% Figure path
param.figdirpath = ['./figs/'];
if ~exist(param.figdirpath) 
    mkdir(param.figdirpath)
end
param.figpath = [param.figdirpath,param.CARDID,'/'];
if ~exist(param.figpath) 
    mkdir(param.figpath)
end

%% Setup mode branches
if param.TYPE == 'S' && branch == 0
    param.periods = param.S0periods;
%     param.periods = [15 20 25 35 50 80 100 150];
    param.Gmat_Spath = param.Gmat_S0path;
elseif param.TYPE == 'S' && branch == 1
    param.periods = param.S1periods;
%     param.periods = [5 6 8 10];
elseif param.TYPE == 'T' && branch == 0
    param.periods = param.T0periods;
%     param.periods = [5 6];
elseif param.TYPE == 'T' && branch == 1 && ~isempty(param.T1periods)
    param.periods = param.T1periods;
%     param.periods = [7];
else
    error('no measurements for this branch!')
end