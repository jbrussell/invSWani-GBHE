function [ data_bs ] = bootstrap_dataset(data, nbs)
% Generate dataset by randomly sampling from a Gaussian distribution with
% mean equal to the observed value and standard deviation equal to the
% measurement uncertainty.
%
%

% Are errors 95% confidence? If so, divide by 2 to get standard deviation.
is_95 = 1;

data_bs(1) = data;
RL_flds = fields(data);
for ibs = 2:nbs
    data_bs(ibs) = data;
    for ifld = 1:length(RL_flds)
        RL = RL_flds{ifld};
        if is_95
            data_bs(ibs).(RL).phi2 = normrnd(data.(RL).phi2,data.(RL).err_phi2/2);
            data_bs(ibs).(RL).phi4 = normrnd(data.(RL).phi4,data.(RL).err_phi4/2);
            data_bs(ibs).(RL).A2 = normrnd(data.(RL).A2,data.(RL).err_2A/2);
            data_bs(ibs).(RL).A4 = normrnd(data.(RL).A4,data.(RL).err_4A/2);
        else
            data_bs(ibs).(RL).phi2 = normrnd(data.(RL).phi2,data.(RL).err_phi2);
            data_bs(ibs).(RL).phi4 = normrnd(data.(RL).phi4,data.(RL).err_phi4);
            data_bs(ibs).(RL).A2 = normrnd(data.(RL).A2,data.(RL).err_2A);
            data_bs(ibs).(RL).A4 = normrnd(data.(RL).A4,data.(RL).err_4A);
        end
    end
end

