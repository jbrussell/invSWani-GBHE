function [ model ] = run_GBHinv( data, model, kernels, par, varargin )
% Invert 2-theta Rayleigh and Love wave phase velocities for G, B, H and 
% following Russell et al., JGR (2019), based off Montagner & Nataf (1986).
%
% If this function called as part of a bootstrap run, chi2 will be 
% calculated with respect to the full dataset. 
%
% 11/19
% jbrussell.github
global data_fds
if ~isempty(varargin)
    ibs = varargin{1};
else
    ibs = [];
end

% Invert cosine 2-theta
[m_c, model] = inv_GBH(data,model,kernels,'c',par,ibs);

% Invert sine 2-theta
[m_s, model] = inv_GBH(data,model,kernels,'s',par,ibs);

% Loop through components and Calculate Strengths and Directions
% FORWARD:
flds = {'R0','R1','L0','L1','R0_fds','R1_fds','L0_fds','L1_fds'};
for ifld = 1:length(flds)
    fld = flds{ifld};
    if isempty(model.(fld).c2)
        model.(fld).A2 = [];
        model.(fld).phi2 = [];
    else
        model.(fld).A2 = sqrt(model.(fld).c2.^2 + model.(fld).s2.^2);
        model.(fld).phi2 = 0.5*atan2d(model.(fld).s2, model.(fld).c2);
        model.(fld).phi2 = find_nearest_azi( model.(fld).phi2, par.FSD, fld )';
    end
end

% MODEL:
flds = {'G','B','H'};
for ifld = 1:length(flds)
    fld = flds{ifld};
    model.(fld).([fld,'_mag']) = sqrt(m_c.(fld).^2 + m_s.(fld).^2);
    model.(fld).([fld,'_dir']) = 0.5*atan2d(m_s.(fld), m_c.(fld));
    model.(fld).([fld,'_dir']) = find_nearest_azi( model.(fld).([fld,'_dir']), par.FSD, fld )';
end

% Calculate chi2
model = model_chi2(model,data,{'R0','R1','L0','L1'});
% % Calculate chi2 (with respect to full dataset)
% model = model_chi2(model,data_fds,{'R0_fds','R1_fds','L0_fds','L1_fds'});

% Combine mode branches into single field
flds = fields(model.R0);
[~,I_sort] = sort([model.R0.periods; model.R1.periods]);
[~,I_sort_fds] = sort([model.R0_fds.periods; model.R1_fds.periods]);
for ifld = 1:length(flds)
    fld = flds{ifld};
    temp = [model.R0.(fld); model.R1.(fld)];
    model.R.(fld) = temp(I_sort);
    if isfield(model.R0_fds,fld)
        temp = [model.R0_fds.(fld); model.R1_fds.(fld)];
        model.R_fds.(fld) = temp(I_sort_fds);
    end
end
model.R.mode_br = data.rayl.mode_br_ani;
model.R_fds.mode_br = data_fds.rayl.mode_br_ani;

flds = fields(model.L0);
[~,I_sort] = sort([model.L0.periods; model.L1.periods]);
[~,I_sort_fds] = sort([model.L0_fds.periods; model.L1_fds.periods]);
for ifld = 1:length(flds)
    fld = flds{ifld};
    temp = [model.L0.(fld); model.L1.(fld)];
    model.L.(fld) = temp(I_sort);
    if isfield(model.L0_fds,fld)
        temp = [model.L0_fds.(fld); model.L1_fds.(fld)];
        model.L_fds.(fld) = temp(I_sort_fds);
    end
end
model.L.mode_br = data.love.mode_br_ani;
model.L_fds.mode_br = data_fds.love.mode_br_ani;

end

function [m, model] = inv_GBH(data,model,kernels,csstr,par,ibs)
    global kernels_fds data_fds
    % Weighted least squares inversion for anisotropic parameters G, B, H
    
    % Initialize some useful things...
    I_R0 = data.rayl.mode_br_ani == 0;
    I_R1 = data.rayl.mode_br_ani == 1;
    I_L0 = data.love.mode_br_ani == 0;
    I_L1 = data.love.mode_br_ani == 1;
    cs2 = [csstr,'2'];
    err_cs2 = ['err_',csstr,'2'];
    G = kernels.G;
    B = kernels.B;
    H = kernels.H;
    nlayerG=G.nlayer; nlayerB=B.nlayer; nlayerH=H.nlayer;
    
    % Setup damping to zero below specified depth: par.dep_zero_damp
    GBH_EYE = eye(G.nlayer+B.nlayer+H.nlayer);
    I_DAMP = find( kernels.G.z >= par.dep_zero_damp);
    G_EYE_DAMP = [GBH_EYE(I_DAMP,:)];
    nlayer_deepDAMP0 = size(G_EYE_DAMP,1);
    % Damp linearly from 0 to par.G_DAMP
    deepDAMP0_vec = linspace(par.G_DAMP*.01,par.G_DAMP,nlayer_deepDAMP0)';
    
    % Damping matrix
    H00 = eye(nlayerG+nlayerB+nlayerH);
    % Smoothing matrix (second derivative)
    [F00_G, F00_B, F00_H] = build_smooth(nlayerG);
    % Flatness matrix (first derivative)
    [J00_G, J00_B, J00_H] = build_flatness(nlayerG);
    
    % Break smoothing constraints by layers?
    if par.is_brk == 1
        [F00_G] = break_constraint(F00_G, G.z, par.z_brks);
        [F00_B] = break_constraint(F00_B, B.z, par.z_brks);
        [F00_H] = break_constraint(F00_H, H.z, par.z_brks);
        [J00_G] = break_constraint(J00_G, G.z, par.z_brks);
        [J00_B] = break_constraint(J00_B, B.z, par.z_brks);
        [J00_H] = break_constraint(J00_H, H.z, par.z_brks);
    end
    
    % Data vector
    d_GBH = [data.rayl.(cs2)(I_R0)'; % R0
             data.rayl.(cs2)(I_R1)'; % R1
             data.love.(cs2)(I_L0)'; % L0
             data.love.(cs2)(I_L1)'];% L1
    
    % Build G matrix
    GGGBBBHHH = [G.K_R0, B.K_R0, H.K_R0; % R0
                 G.K_R1, B.K_R1, H.K_R1; % R1
                 G.K_L0, zeros(size(G.K_L0)), zeros(size(G.K_L0)); % L0
                 G.K_L1, zeros(size(G.K_L1)), zeros(size(G.K_L1));]; % L1
    
    % Data weights
    Wd = diag([1./data.rayl.(err_cs2)(I_R0)'; % R0
               1./data.rayl.(err_cs2)(I_R1)'; % R1
               1./data.love.(err_cs2)(I_L0)'; % L0
               1./data.love.(err_cs2)(I_L1)']).^2; % L1
    % Downweight large weights 
    Wd(Wd>2*std(diag(Wd))) = 2*std(diag(Wd));
    
    % Constraint equations
    h = [zeros(nlayer_deepDAMP0,1); % G damping to zero below specified depth
         zeros(nlayerB,1); % B/G ratio
         zeros(nlayerH,1); % H/G ratio
         zeros(size(H00,1),1); % G,B,H Norm damping
         zeros(size(F00_G,1),1); % G 2nd derivative smoothing
         zeros(size(F00_B,1),1); % B 2nd derivative smoothing
         zeros(size(F00_H,1),1); % H 2nd derivative smoothing
         zeros(size(J00_G,1),1); % G 1st derivative smoothing
         zeros(size(J00_B,1),1); % B 1st derivative smoothing
         zeros(size(J00_H,1),1); % H 1st derivative smoothing
         ];
     
    H_mat = [G_EYE_DAMP;  % G damping to zero below specified depth
        -par.BGratio*eye(nlayerG), eye(nlayerB), zeros(nlayerH); % B/G ratio
        -par.HGratio*eye(nlayerG), zeros(nlayerB), eye(nlayerH); % H/G ratio
         H00; % G,B,H Norm damping
         F00_G; % G 2nd derivative smoothing
         F00_B; % B 2nd derivative smoothing
         F00_H; % H 2nd derivative smoothing
         J00_G; % G 1st derivative smoothing
         J00_B; % B 1st derivative smoothing
         J00_H; % H 1st derivative smoothing
         ];
   
    % Constraint equation weights
    Wh = diag([ones(nlayer_deepDAMP0,1).*deepDAMP0_vec; % G damping to zero below specified depth
               ones(nlayerG,1)*par.epsilonBG; % B/G ratio
               ones(nlayerG,1)*par.epsilonHG; % H/G ratio
               ones(size(H00,1),1)*par.alphH_GBH; % G,B,H Norm damping
               ones(size(F00_G,1),1)*par.alphF_GBH; % G 2nd derivative smoothing
               ones(size(F00_B,1),1)*par.alphF_GBH; % B 2nd derivative smoothing
               ones(size(F00_H,1),1)*par.alphF_GBH; % H 2nd derivative smoothing
               ones(size(J00_G,1),1)*par.alphJ_GBH; % G 1st derivative smoothing
               ones(size(J00_B,1),1)*par.alphJ_GBH; % B 1st derivative smoothing
               ones(size(J00_H,1),1)*par.alphJ_GBH; % H 1st derivative smoothing
               ]).^2;
           
    % Setup augmented matrix
    f = [par.eps_d*Wd.^(1/2)*d_GBH; % data vector
         par.eps_f*Wh.^(1/2)*h]; % constraint vector
    
    F = [par.eps_d*Wd.^(1/2)*GGGBBBHHH; % G matrix for data
         par.eps_f*Wh.^(1/2)*H_mat]; % constraint equations
            
    % Least squares solution
    m_GBH = (F'*F)\F'*f;
    
    % Calculate effective number of model parameters, data, and degrees of freedom
    GGGBBBHHHinv = (F'*F)\GGGBBBHHH'*Wd;
    % Model resolution
    R = GGGBBBHHHinv * GGGBBBHHH;
%     R2 = (GGGBBBHHH'*Wd*GGGBBBHHH + H_mat'*Wh*H_mat)\GGGBBBHHH'*Wd*GGGBBBHHH;
    % Data resolution
    N = GGGBBBHHH * GGGBBBHHHinv;
%     N2 = GGGBBBHHH*((GGGBBBHHH'*Wd*GGGBBBHHH + H_mat'*Wh*H_mat)\GGGBBBHHH'*Wd);
    % Effective degrees of freedom
    model.chi2.(['df2_eff_',cs2]) = length(d_GBH) - trace(R);
    if  isempty(ibs) || ibs == 1
        model.resolution.(['R_GBH',cs2]) = R;
        model.resolution.(['N_GBH',cs2]) = N;
    end
    R_trace = cumsum(diag(R));
    model.resolution.(['traceR_GBH',cs2]) = [R_trace(1:length(G.z)) G.z];

    % Index G, B, and H
    m_G = m_GBH(1:nlayerG);
    m_B = m_GBH(nlayerG+1:nlayerG+nlayerB);
    m_H = m_GBH(nlayerG+nlayerB+1:nlayerG+nlayerB+nlayerH);
    % display(m_B./m_G);
    % display(m_H./m_G);
    
    % Forward calculate observations
    model.R0.(cs2) = G.K_R0*m_G + B.K_R0*m_B + H.K_R0*m_H;
    model.R1.(cs2) = G.K_R1*m_G + B.K_R1*m_B + H.K_R1*m_H;
    model.L0.(cs2) = G.K_L0*m_G;
    if ~isempty(G.K_L1)
        model.L1.(cs2) = G.K_L1*m_G;
    else
        model.L1.(cs2) = [];
    end
    % (on full dataset)
    model.R0_fds.(cs2) = kernels_fds.G.K_R0*m_G + kernels_fds.B.K_R0*m_B + kernels_fds.H.K_R0*m_H;
    model.R1_fds.(cs2) = kernels_fds.G.K_R1*m_G + kernels_fds.B.K_R1*m_B + kernels_fds.H.K_R1*m_H;
    model.L0_fds.(cs2) = kernels_fds.G.K_L0*m_G;
    if ~isempty(kernels_fds.G.K_L1)
        model.L1_fds.(cs2) = kernels_fds.G.K_L1*m_G;
    else
        model.L1_fds.(cs2) = [];
    end
    
    % Calculate residuals
    model.R0.(['d',cs2]) = data.rayl.(cs2)(I_R0)' - model.R0.(cs2);
    model.R1.(['d',cs2]) = data.rayl.(cs2)(I_R1)' - model.R1.(cs2);
    model.L0.(['d',cs2]) = data.love.(cs2)(I_L0)' - model.L0.(cs2);
    if ~isempty(G.K_L1)
        model.L1.(['d',cs2]) = data.love.(cs2)(I_L1)' - model.L1.(cs2);
    else
        model.L1.(['d',cs2]) = [];
    end
    
    % (on full dataset)
    % First add "periods" field to fds
    I.R0_fds = data_fds.rayl.mode_br_ani == 0;
    I.R1_fds = data_fds.rayl.mode_br_ani == 1;
    I.L0_fds = data_fds.love.mode_br_ani == 0;
    I.L1_fds = data_fds.love.mode_br_ani == 1;
    flds = {'R0','R1','L0','L1'};
    for ifld = 1:length(flds)
        if strcmp(flds{ifld}(1),'R')
            type = 'rayl';
        else
            type = 'love';
        end
        model.([flds{ifld},'_fds']).periods = data_fds.(type).periods_ani(I.([flds{ifld},'_fds']))';
    end
    model.R0_fds.(['d',cs2]) = data_fds.rayl.(cs2)(I.R0_fds)' - model.R0_fds.(cs2);
    model.R1_fds.(['d',cs2]) = data_fds.rayl.(cs2)(I.R1_fds)' - model.R1_fds.(cs2);
    model.L0_fds.(['d',cs2]) = data_fds.love.(cs2)(I.L0_fds)' - model.L0_fds.(cs2);
    if ~isempty(kernels_fds.G.K_L1)
        model.L1_fds.(['d',cs2]) = data_fds.love.(cs2)(I.L1_fds)' - model.L1_fds.(cs2);
    else
        model.L1_fds.(['d',cs2]) = [];
    end
    
    
    model.G.(['G',csstr]) = m_G;
    model.B.(['B',csstr]) = m_B;
    model.H.(['H',csstr]) = m_H;
    
    m.G = m_G;
    m.B = m_B;
    m.H = m_H;
end

function [F00_G, F00_B, F00_H] = build_smooth(nlayer)
    % smoothing matrix
    F00 = 2*eye(1*nlayer);
    Fup = -1*[zeros(1*nlayer-1,1) eye(1*nlayer-1) ;zeros(1,1*nlayer) ];
    Fdown = -1*[zeros(1,1*nlayer); eye(1*nlayer-1) zeros(1*nlayer-1,1) ];
    F00 = F00+Fup+Fdown;
    F00(1,:) = 0; F00(end,:) = 0;
    F00_G = [F00, zeros(nlayer), zeros(nlayer)];
    F00_B = [zeros(nlayer), F00, zeros(nlayer)];
    F00_H = [zeros(nlayer), zeros(nlayer), F00];
end

function [J00_G, J00_B, J00_H] = build_flatness(nlayer)
    % Flatness matrix (first derviative)
    J00 = 1*eye(nlayer);
    Jdown = -1*[zeros(1,nlayer); eye(nlayer-1) zeros(nlayer-1,1) ];
    J00 = J00+Jdown;
    J00(1,:) = 0; J00(end,:) = 0;
    J00_G = [J00, zeros(nlayer), zeros(nlayer)];
    J00_B = [zeros(nlayer), J00, zeros(nlayer)];
    J00_H = [zeros(nlayer), zeros(nlayer), J00];
end

function [ofastdirs] = find_nearest_azi(nfastdirs,azi,fld)
    % Find fast direction nearest to azi
    for idep = 1:length(nfastdirs)    
        fastdir_vec(1) = nfastdirs(idep);
        fastdir_vec(2) = nfastdirs(idep)+180;
        fastdir_vec(3) = nfastdirs(idep)-180;
        if strcmp(fld,'H') == 1
            [~, I] = min(abs(fastdir_vec-azi-90));
        else
            [~, I] = min(abs(fastdir_vec-azi));
        end
        fastdir_hold(idep) = fastdir_vec(I);
    end
    ind = find(fastdir_hold<0);
    fastdir_hold(ind) = fastdir_hold(ind)+180;
    ofastdirs = fastdir_hold;
end

function [model] = model_chi2(model,data,flds)
    % Initialize some useful things...
    I.(flds{1}) = data.rayl.mode_br_ani == 0;
    I.(flds{2}) = data.rayl.mode_br_ani == 1;
    I.(flds{3}) = data.love.mode_br_ani == 0;
    I.(flds{4}) = data.love.mode_br_ani == 1;
    model.chi2.X2_2 = 0;
    model.chi2.df2 = 0;
    for ifld = 1:length(flds)
        fld = flds{ifld};
        if sum(I.(fld))==0
            continue
        end
        if strcmp(fld(1),'R')
            d = data.rayl;
        elseif strcmp(fld(1),'L')
            d = data.love;
        end
        model.chi2.(fld).chi2_2 = sum(model.(fld).dc2.^2./d.err_c2(I.(fld)).^2') + sum(model.(fld).ds2.^2./d.err_s2(I.(fld)).^2');
        model.chi2.X2_2 = model.chi2.X2_2 + model.chi2.(fld).chi2_2;
        model.chi2.df2 = model.chi2.df2 + length(d.err_c2(I.(fld))) + length(d.err_s2(I.(fld)));
    end
    model.chi2.df2_eff = model.chi2.df2_eff_c2 + model.chi2.df2_eff_s2;
%     model.chi2.X2red_2 = model.chi2.X2_2 / model.chi2.df2;
    model.chi2.X2red_2 = model.chi2.X2_2 / model.chi2.df2_eff;
end

function [X00] = break_constraint(X00, z, z_brks)
% Break constraint equation at specified depths
    for ibrk = 1:length(z_brks)
        if z_brks(ibrk) > max(z)
            continue
        end
        [~,I_brk] = min(abs(z-z_brks(ibrk)));
        if length(find(X00(3,:)~=0))==3 % second derivative
            X00(I_brk-1:I_brk+1,:) = 0;
        end
        if length(find(X00(3,:)~=0))==2 % first derivative
            X00(I_brk,:) = 0;
        end
    end
end
