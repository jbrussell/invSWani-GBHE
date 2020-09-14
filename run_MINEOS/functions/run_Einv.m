function [ model ] = run_Einv( data, model, kernels, par, varargin )
% Do the E inversion following Russell et al., JGR (2019), based off
% Montagner & Nataf (1986).
%
% 11/19
% jbrussell.github
global data_fds
if ~isempty(varargin)
    ibs = varargin{1};
else
    ibs = [];
end

% Invert cosine 4-theta
[m_c, model] = inv_E(data,model,kernels,'c',par,ibs);

% Invert sine 4-theta
[m_s, model] = inv_E(data,model,kernels,'s',par,ibs);

% Calculate Strength and Direction
% FORWARD:
flds = {'L0','L1','L0_fds','L1_fds'};
for ifld = 1:length(flds)
    fld = flds{ifld};
    if isempty(model.(fld).c4)
        model.(fld).A4 = [];
        model.(fld).phi4 = [];
    else
        model.(fld).A4 = sqrt(model.(fld).c4.^2 + model.(fld).s4.^2);
        model.(fld).phi4 = 0.25*atan2d(model.(fld).s4, model.(fld).c4);
        model.(fld).phi4 = find_nearest_azi( model.(fld).phi4, par.FSD )';
    end
end

% MODEL:
model.E.E_mag = sqrt(m_c.E.^2 + m_s.E.^2);
model.E.E_dir = 0.25*atan2d(-m_s.E, -m_c.E);
model.E.E_dir = find_nearest_azi( model.E.E_dir, par.FSD )';

% Calculate chi2
model = model_chi2(model,data,{'L0','L1'});
% model = model_chi2(model,data_fds,{'L0_fds','L1_fds'});

% Combine mode branches into single field
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

function [m, model] = inv_E(data,model,kernels,csstr,par,ibs)
    global kernels_fds data_fds
    % Weighted least squares inversion for anisotropic parameter E. Invert
    % following Menke (2012), Environmental Data Analysis; Ch. 5.5; equation 5.17

    % Initialize some useful things...
    I_R0 = data.rayl.mode_br_ani == 0;
    I_R1 = data.rayl.mode_br_ani == 1;
    I_L0 = data.love.mode_br_ani == 0;
    I_L1 = data.love.mode_br_ani == 1;
    cs4 = [csstr,'4'];
    err_cs4 = ['err_',csstr,'4'];
    E = kernels.E;
    nlayerE=E.nlayer;
    
    % Damping matrix
    H0 = eye(nlayerE);
    % Smoothing matrix
    [F0] = build_smooth(nlayerE);
    % Flatness matrix
    [J0] = build_flatness(nlayerE);
    
    % Break smoothing constraints by layers?
    if par.is_brk == 1
        [F0] = break_constraint(F0, E.z, par.z_brks);
        [J0] = break_constraint(J0, E.z, par.z_brks);
    end
            
    % Data vector
    d_E = [data.love.(cs4)(I_L0)'; % L0
           data.love.(cs4)(I_L1)'];% L1
    
    % Build G matrix
    EEE = [E.K_L0; % L0
           E.K_L1; % L1
           ];
       
    % Data weights
    Wd = diag([1./data.love.(err_cs4)(I_L0)'; % L0
               1./data.love.(err_cs4)(I_L1)']).^2; % L1
    % Downweight large weights 
    Wd(Wd>2*std(diag(Wd))) = 2*std(diag(Wd));
           
   % Constraint equations
    h = [zeros(size(H0,1),1); % E Norm damping
         zeros(size(J0,1),1); % E 1st derivative smoothing
         zeros(size(F0,1),1); % E 2nd derivative smoothing
         ];
     
    H_mat = [H0; % E Norm damping
             J0; % E First derivative smoothing
             F0; % E Second derivative smoothing
             ];
   
    % Constraint equation weights
    Wh = diag([ones(size(H0,1),1)*par.alphH_E; % E Norm damping
               ones(size(J0,1),1)*par.alphJ_E; % E 1st derivative smoothing
               ones(size(F0,1),1)*par.alphF_E; % E 2nd derivative smoothing
               ]).^2;
           
    % Setup augmented matrix
    f = [par.eps_d*Wd.^(1/2)*d_E; % data vector
         par.eps_f*Wh.^(1/2)*h]; % constraint vector
    
    F = [par.eps_d*Wd.^(1/2)*EEE; % G matrix for data
         par.eps_f*Wh.^(1/2)*H_mat]; % constraint equations
            
    % Least squares solution
    m_E = (F'*F)\F'*f;
    
    % Calculate effective number of model parameters, data, and degrees of freedom
    EEEinv = (F'*F)\EEE'*Wd;
    % Model resolution
    R = EEEinv * EEE;
    % Data resolution
    N = EEE * EEEinv;
    % Effective degrees of freedom
    model.chi2.(['df4_eff_',cs4]) = length(d_E) - trace(R);
    if isempty(ibs) || ibs == 1
        model.resolution.(['R_E',cs4]) = R;
        model.resolution.(['N_E',cs4]) = N;
    end
    R_trace = cumsum(diag(R));
    model.resolution.(['traceR_E',cs4]) = [R_trace(1:length(E.z)) E.z];
    
    % Forward calculate observations
    model.L0.(cs4) = E.K_L0*m_E;
    if ~isempty(E.K_L1)
        model.L1.(cs4) = E.K_L1*m_E;
    else
        model.L1.(cs4) = [];
    end
    
    % (on full dataset)
    model.L0_fds.(cs4) = kernels_fds.E.K_L0*m_E;
    if ~isempty(kernels_fds.E.K_L1)
        model.L1_fds.(cs4) = kernels_fds.E.K_L1*m_E;
    else
        model.L1_fds.(cs4) = [];
    end
    
    % Calculate residuals
    model.L0.(['d',cs4]) = data.love.(cs4)(I_L0)' - model.L0.(cs4);
    if ~isempty(E.K_L1)
        model.L1.(['d',cs4]) = data.love.(cs4)(I_L1)' - model.L1.(cs4);
    else
        model.L1.(['d',cs4]) = [];
    end
    
    % (on full dataset) 
    % First add "periods" field to fds
    I.L0_fds = data_fds.love.mode_br_ani == 0;
    I.L1_fds = data_fds.love.mode_br_ani == 1;
    flds = {'L0','L1'};
    for ifld = 1:length(flds)
        model.([flds{ifld},'_fds']).periods = data_fds.love.periods_ani(I.([flds{ifld},'_fds']))';
    end    
    model.L0_fds.(['d',cs4]) = data_fds.love.(cs4)(I.L0_fds)' - model.L0_fds.(cs4);
    if ~isempty(kernels_fds.E.K_L1)
        model.L1_fds.(['d',cs4]) = data_fds.love.(cs4)(I.L1_fds)' - model.L1_fds.(cs4);
    else
        model.L1_fds.(['d',cs4]) = [];
    end
    
    model.E.(['E',csstr]) = m_E;
    
    m.E = m_E;
end

function [F00] = build_smooth(nlayer)
    % smoothing matrix
    F00 = 2*eye(1*nlayer);
    Fup = -1*[zeros(1*nlayer-1,1) eye(1*nlayer-1) ;zeros(1,1*nlayer) ];
    Fdown = -1*[zeros(1,1*nlayer); eye(1*nlayer-1) zeros(1*nlayer-1,1) ];
    F00 = F00+Fup+Fdown;
    F00(1,:) = 0; F00(end,:) = 0;
    F00(1,1) = 1; F00(1,2) = -1;
    F00(end,end-1) = 1; F00(end,end) = -1;
end

function [J0] = build_flatness(nlayer)
    % Flatness matrix (first derviative)
    J0 = 1*eye(nlayer);
    Jdown = -1*[zeros(1,nlayer); eye(nlayer-1) zeros(nlayer-1,1) ];
    J0 = J0+Jdown;
    J0(1,:) = 0; J0(end,:) = 0;
end

function [ofastdirs] = find_nearest_azi(nfastdirs,azi)
    % Find fast direction nearest to azi
    ind = find(nfastdirs<0);
    nfastdirs(ind) = nfastdirs(ind)+180;
    for idep = 1:length(nfastdirs)    
        fastdir_vec(1) = nfastdirs(idep);
        fastdir_vec(2) = nfastdirs(idep)+90;
        fastdir_vec(3) = nfastdirs(idep)+180;
        fastdir_vec(4) = nfastdirs(idep)+270;
        fastdir_vec(5) = nfastdirs(idep)-90;
        fastdir_vec(6) = nfastdirs(idep)-180;
        fastdir_vec(7) = nfastdirs(idep)-270;
        [~, I] = min(abs(fastdir_vec-azi-45));
        fastdir_hold(idep) = fastdir_vec(I);       
    end
    ofastdirs = fastdir_hold;
end

function [model] = model_chi2(model,data,flds)
    % Initialize some useful things...
%     I.R0 = data.rayl.mode_br_ani == 0;
%     I.R1 = data.rayl.mode_br_ani == 1;
    I.(flds{1}) = data.love.mode_br_ani == 0;
    I.(flds{2}) = data.love.mode_br_ani == 1;
    model.chi2.X2_4 = 0;
    model.chi2.df4 = 0;
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
        model.chi2.(fld).chi2_4 = sum(model.(fld).dc4.^2./d.err_c4(I.(fld)).^2') + sum(model.(fld).ds4.^2./d.err_s4(I.(fld)).^2');
        model.chi2.X2_4 = model.chi2.X2_4 + model.chi2.(fld).chi2_4;
        model.chi2.df4 = model.chi2.df4 + length(d.err_c4(I.(fld))) + length(d.err_s4(I.(fld)));
    end
    model.chi2.df4_eff = model.chi2.df4_eff_c4 + model.chi2.df4_eff_s4;
%     model.chi2.X2red_4 = model.chi2.X2_4 / model.chi2.df4;
    model.chi2.X2red_4 = model.chi2.X2_4 / model.chi2.df4;
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