function [ fwd ] = forward_GBH( fwd,kernels,data )
% Forward calculate R0, R1, and L phase velocity data from input 
% (G_mag,G_dir), (B_mag,B_dir), (H_mag,H_dir) models
%
% INPUT:
%   - fwd.G.G_mag
%   - fwd.G.G_dir
%   - fwd.B.B_mag
%   - fwd.B.B_dir
%   - fwd.H.H_mag
%   - fwd.H.H_dir
%   - kernels structure output from load_kernels.m
%
% JBR: 12/19

% Convert (Mag,Dir) to (sin,cos)
flds_GBH = {'G','B','H'};
for ifld = 1:length(flds_GBH)
    GBH = flds_GBH{ifld};
    fwd.(GBH).([GBH,'c']) = fwd.(GBH).([GBH,'_mag']).*cosd(2*fwd.(GBH).([GBH,'_dir']));
    fwd.(GBH).([GBH,'s']) = fwd.(GBH).([GBH,'_mag']).*sind(2*fwd.(GBH).([GBH,'_dir']));
end

% Forward calculate phase velocities using kernels
flds_RL = {'R0','R1','L0','L1'};
for ifld = 1:length(flds_RL)
    RL = flds_RL{ifld};
    if isempty(kernels.G.(['K_',RL]))
        fwd.(RL).c2=[];  fwd.(RL).s2=[];
        continue
    end    
    % Loop through kernels and sum
    fwd.(RL).c2=0;  fwd.(RL).s2=0;
    for jfld = 1:length(flds_GBH)
        GBH = flds_GBH{jfld};
        GBHc=[GBH,'c'];  GBHs=[GBH,'s'];
        if ~isfield(kernels.(GBH),['K_',RL])
            continue
        end
        fwd.(RL).c2 = fwd.(RL).c2 + kernels.(GBH).(['K_',RL])*fwd.(GBH).(GBHc);
        fwd.(RL).s2 = fwd.(RL).s2 + kernels.(GBH).(['K_',RL])*fwd.(GBH).(GBHs);
    end
end

% Convert (sin,cos) back to (A,phi)
for ifld = 1:length(flds_RL)
    RL = flds_RL{ifld};
    if isempty(fwd.(RL).c2)
        fwd.(RL).A2=[];  fwd.(RL).phi2=[];
        continue
    end
    fwd.(RL).A2 = sqrt(fwd.(RL).c2.^2 + fwd.(RL).s2.^2);
    fwd.(RL).phi2 = 0.5*atan2d(fwd.(RL).s2, fwd.(RL).c2);
    fwd.(RL).phi2(fwd.(RL).phi2<0) = fwd.(RL).phi2(fwd.(RL).phi2<0) + 180;
end

% Calculate chi2 misfit
I.R0 = data.rayl.mode_br_ani == 0;
I.R1 = data.rayl.mode_br_ani == 1;
I.L0 = data.love.mode_br_ani == 0;
I.L1 = data.love.mode_br_ani == 1;
fwd.R0.dc2 = data.rayl.c2(I.R0)' - fwd.R0.c2;
fwd.R1.dc2 = data.rayl.c2(I.R1)' - fwd.R1.c2;
fwd.R0.ds2 = data.rayl.s2(I.R0)' - fwd.R0.s2;
fwd.R1.ds2 = data.rayl.s2(I.R1)' - fwd.R1.s2;
fwd.R0.dA2 = data.rayl.A2(I.R0)' - fwd.R0.A2;
fwd.R1.dA2 = data.rayl.A2(I.R1)' - fwd.R1.A2;
fwd.R0.dphi2 = data.rayl.phi2(I.R0)' - fwd.R0.phi2;
fwd.R1.dphi2 = data.rayl.phi2(I.R1)' - fwd.R1.phi2;
fwd.L0.dc2 = data.love.c2(I.L0)' - fwd.L0.c2;
fwd.L1.dc2 = data.love.c2(I.L1)' - fwd.L1.c2;
fwd.L0.ds2 = data.love.s2(I.L0)' - fwd.L0.s2;
fwd.L1.ds2 = data.love.s2(I.L1)' - fwd.L1.s2;
fwd = model_chi2(fwd,data,{'R0','R1','L0','L1'});


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
%         model.chi2.(fld).chi2_2 = sum(model.(fld).dA2.^2./d.err_2A(I.(fld)).^2') + sum(model.(fld).dphi2.^2./d.err_phi2(I.(fld)).^2');
        model.chi2.X2_2 = model.chi2.X2_2 + model.chi2.(fld).chi2_2;
        model.chi2.df2 = model.chi2.df2 + length(d.err_c2(I.(fld))) + length(d.err_s2(I.(fld)));
    end
%     model.chi2.df2_eff = model.chi2.df2_eff_c2 + model.chi2.df2_eff_s2;
    model.chi2.X2red_2 = model.chi2.X2_2 / (model.chi2.df2-10);
%     model.chi2.X2red_2 = model.chi2.X2_2 / model.chi2.df2_eff;
end

