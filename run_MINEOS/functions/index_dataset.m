function [ data_bs ] = index_dataset( data, nbs )
% Index dataset for balanced bootstrap resampling
%
% JBR: 11/19

% Get number of measurements
ndata_R0 = sum(data.rayl.mode_br_ani==0);
ndata_R1 = sum(data.rayl.mode_br_ani==1);
ndata_R = ndata_R0 + ndata_R1;
ndata_L0 = sum(data.love.mode_br_ani==0);
ndata_L1 = sum(data.love.mode_br_ani==1);
ndata_L = ndata_L0 + ndata_L1;

% Setup indeces for balanced resampling
[ indxs_R ] = balanced_resampling( ndata_R,nbs-1 );
indxs_R = sort(indxs_R,1);
indxs_R = [[1:ndata_R]' , indxs_R];

[ indxs_L ] = balanced_resampling( ndata_L,nbs-1 );
indxs_L = sort(indxs_L,1);
indxs_L = [[1:ndata_L]' , indxs_L];

% [ indxs_RL ] = balanced_resampling( ndata_R + ndata_L,nbs-1 );
% indxs_RL = sort(indxs_RL,1);
% indxs_RL = [[1:ndata_R + ndata_L]' , indxs_RL];

% Index Rayleigh and Love measurements
R_flds = fields(data.rayl);
L_flds = fields(data.love);
for ibs = 1:size(indxs_R,2)
    for ifld = 1:length(R_flds)
        data_bs(ibs).rayl.(R_flds{ifld}) = data.rayl.(R_flds{ifld})(indxs_R(:,ibs));
    end
    for ifld = 1:length(L_flds)
        data_bs(ibs).love.(L_flds{ifld}) = data.love.(L_flds{ifld})(indxs_L(:,ibs));
    end
    data_bs(ibs).indx_R = indxs_R(:,ibs);
    data_bs(ibs).indx_L = indxs_L(:,ibs);
end

