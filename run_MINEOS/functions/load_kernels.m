function [ kernels ] = load_kernels( param, par, data, model )
% Load in sensitivity kernels, organize them by frequency, and scale them
% according to Russell et al., JGR (2019)
%
% Kernel dimensions: (F x Z)
% F - number of frequencies, increasing from top to bottom
% Z - number of layers, incresing depth to the right.
% 
%
% 11/19
% jbrussell.github

Vsv = model.card_midpts.Vsv;  Vsh = model.card_midpts.Vsh;
Vpv = model.card_midpts.Vpv;  Vph = model.card_midpts.Vph;
rho = model.card_midpts.rho;  eta = model.card_midpts.eta;
rho_E = model.card_midpts.rho_E; Vsh_E = model.card_midpts.Vsh_E;
Gmat_S0path = param.Gmat_S0path;
Gmat_S1path = param.Gmat_S1path;
Gmat_T0path = param.Gmat_T0path;
if ~isempty(data.love.periods_ani(data.love.mode_br_ani == 1))
    isT1 = 1;
    Gmat_T1path = param.Gmat_T1path;
else
    isT1 = 0;
end

%% Load kernels
% Load G matrices
T0_Gmat = load(Gmat_T0path);
S0_Gmat = load(Gmat_S0path);
S1_Gmat = load(Gmat_S1path);
if isT1
    T1_Gmat = load(Gmat_T1path);
else
    T1_Gmat = [];
end

% Ensure only necessary frequencies are used
[ S0_Gmat, S1_Gmat, T0_Gmat, T1_Gmat ] = clean_kernels( S0_Gmat, S1_Gmat, T0_Gmat, T1_Gmat, data );

% Kernel depth indices
I = find(model.card.vpv > par.Vpv_cutoff);
qeqeqe = length(model.card.vpv)-I(end);
G.Idep = find(S0_Gmat.Gmatrix.rad >= S0_Gmat.Gmatrix.rad(1)-par.model_depth_G*1000);
G.Idep = G.Idep(qeqeqe+1:end);
B.Idep = find(S0_Gmat.Gmatrix.rad >= S0_Gmat.Gmatrix.rad(1)-par.model_depth_B*1000);
B.Idep = B.Idep(qeqeqe+1:end);
H.Idep = find(S0_Gmat.Gmatrix.rad >= S0_Gmat.Gmatrix.rad(1)-par.model_depth_H*1000);
H.Idep = H.Idep(qeqeqe+1:end);
E.Idep = find(T0_Gmat.Gmatrix.rad >= T0_Gmat.Gmatrix.rad(1)-par.model_depth_E*1000);
E.Idep = E.Idep(qeqeqe+1:end);

% Calculate dr terms for multiplying matrix
radG = S0_Gmat.Gmatrix.rad(G.Idep);
G.dr = -[diff(radG)]';
B.dr = G.dr;
H.dr = G.dr;
radE = T0_Gmat.Gmatrix.rad(E.Idep);
E.dr = -[diff(radE)]';

G.rad = radG;  G.z = midpts(S0_Gmat.Gmatrix.rad(1)- G.rad')'/1000;
B.rad = radG;  B.z = midpts(S0_Gmat.Gmatrix.rad(1)- B.rad')'/1000;
H.rad = radG;  H.z = midpts(S0_Gmat.Gmatrix.rad(1)- H.rad')'/1000;
E.rad = radE;  E.z = midpts(T0_Gmat.Gmatrix.rad(1)- E.rad')'/1000;

G.nlayer = length(G.dr);
B.nlayer = length(B.dr);
H.nlayer = length(H.dr);
E.nlayer = length(E.dr);

% Get periods
G.periodsR0 = S0_Gmat.Gmatrix.periods;
G.periodsR1 = S1_Gmat.Gmatrix.periods;
G.periodsL0 = T0_Gmat.Gmatrix.periods;
if isT1
    G.periodsL1 = T1_Gmat.Gmatrix.periods;
else
    G.periodsL1 = [];
end
B.periodsR0 = S0_Gmat.Gmatrix.periods;
B.periodsR1 = S1_Gmat.Gmatrix.periods;
H.periodsR0 = S0_Gmat.Gmatrix.periods;
H.periodsR1 = S1_Gmat.Gmatrix.periods;
E.periodsL0 = T0_Gmat.Gmatrix.periods;
if isT1
    E.periodsL1 = T1_Gmat.Gmatrix.periods;
else
    E.periodsL1 = [];
end


%% Save unscaled kernels (Montagner & Nataf; 1986)
% 2-theta terms
% G: (L) Rayleigh, (-L) Love
G.K_R0raw = midpts(S0_Gmat.Gmatrix.L(:,G.Idep));
G.K_R1raw = midpts(S1_Gmat.Gmatrix.L(:,G.Idep));
G.K_L0raw = midpts(-T0_Gmat.Gmatrix.L(:,G.Idep));
if isT1
    G.K_L1raw = midpts(-T1_Gmat.Gmatrix.L(:,G.Idep));
else
    G.K_L1raw = [];
end
% B: (A) Rayleigh
B.K_R0raw = midpts(S0_Gmat.Gmatrix.A(:,B.Idep));
B.K_R1raw = midpts(S1_Gmat.Gmatrix.A(:,B.Idep));
% H: (F) Rayleigh
H.K_R0raw = midpts(S0_Gmat.Gmatrix.F(:,B.Idep));
H.K_R1raw = midpts(S1_Gmat.Gmatrix.F(:,B.Idep));

% 4-theta terms
% E: (-N) Love, (A) Rayleigh [ignored]
E.K_L0raw = midpts(-T0_Gmat.Gmatrix.N(:,E.Idep));
if isT1
    E.K_L1raw = midpts(-T1_Gmat.Gmatrix.N(:,E.Idep));
else
    E.K_L1raw = [];
end

%% Calculate scaling terms
G.scale_R0 = (model.R0.phv./model.R0.grv) .* rho' .* Vsv'.^2;
G.scale_R1 = (model.R1.phv./model.R1.grv) .* rho' .* Vsv'.^2;
G.scale_L0 = (model.L0.phv./model.L0.grv) .* rho' .* Vsv'.^2;
if isT1
    G.scale_L1 = (model.L1.phv./model.L1.grv) .* rho' .* Vsv'.^2;
else
    G.scale_L1 = [];
end

B.scale_R0 = (model.R0.phv./model.R0.grv) .* rho' .* Vph'.^2;
B.scale_R1 = (model.R1.phv./model.R1.grv) .* rho' .* Vph'.^2;

H.scale_R0 = (model.R0.phv./model.R0.grv) .* rho' .*eta' .* (Vph'.^2 - 2*Vsv'.^2);
H.scale_R1 = (model.R1.phv./model.R1.grv) .* rho' .*eta' .* (Vph'.^2 - 2*Vsv'.^2);

E.scale_L0 = (model.L0.phv./model.L0.grv) .* rho_E' .* Vsh_E'.^2;
if isT1
    E.scale_L1 = (model.L1.phv./model.L1.grv) .* rho_E' .* Vsh_E'.^2;
else
    E.scale_L1 = [];
end


%% Apply scaling
G.K_R0 = G.K_R0raw .* G.scale_R0 .* G.dr;
G.K_R1 = G.K_R1raw .* G.scale_R1 .* G.dr;
G.K_L0 = G.K_L0raw .* G.scale_L0 .* G.dr;
if isT1
    G.K_L1 = G.K_L1raw .* G.scale_L1 .* G.dr;
else
    G.K_L1 = [];
end

B.K_R0 = B.K_R0raw .* B.scale_R0 .* B.dr;
B.K_R1 = B.K_R1raw .* B.scale_R1 .* B.dr;

H.K_R0 = H.K_R0raw .* H.scale_R0 .* H.dr;
H.K_R1 = H.K_R1raw .* H.scale_R1 .* H.dr;

E.K_L0 = E.K_L0raw .* E.scale_L0 .* E.dr;
if isT1
    E.K_L1 = E.K_L1raw .* E.scale_L1 .* E.dr;
else
    E.K_L1 = [];
end

kernels.G = G;
kernels.B = B;
kernels.H = H;
kernels.E = E;
end

