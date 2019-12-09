function [ oS0_Gmat, oS1_Gmat, oT0_Gmat, oT1_Gmat ] = clean_kernels( S0_Gmat, S1_Gmat, T0_Gmat, T1_Gmat, data )
% Organize kernels such that only necessary frequencies are used.

% Initialize
Gmats = {S0_Gmat, S1_Gmat, T0_Gmat, T1_Gmat};
modes = [0, 1, 0, 1];
types = {'rayl','rayl','love','love'};

% Ensure that kernel frequencies match data
for iG = 1:length(Gmats)
    if isempty(Gmats{iG})
        continue
    end
    Gmat = Gmats{iG}.Gmatrix;
    
    % Loop over periods and find those that match data
    periods_data = data.(types{iG}).periods_ani;  
    periods_data = periods_data( data.(types{iG}).mode_br_ani==modes(iG) );
    ii = 0;
    I_save = 0;
    for ip = 1:length(periods_data)
        if ismember(periods_data(ip),Gmat.periods)
            ii = ii + 1;
%             I_save(ii) = ip;
            I_save(ii) = find(Gmat.periods==periods_data(ip));
        end
    end
    
    % Loop over fields and index correct values
    flds = fields(Gmat);
    for ifld = 1:length(flds)
        if ismember(flds{ifld},{'A','C','F','L','N'})
            Gmat.(flds{ifld}) = Gmat.(flds{ifld})(I_save,:);
        end
    end
    Gmat.periods = Gmat.periods(I_save);
    Gmats{iG}.Gmatrix = Gmat;
    
%     if length(Gmat.periods) ~= length(periods_data)
%         error('Number of kernels does not match number of data')
%     end
end

oS0_Gmat = Gmats{1};
oS1_Gmat = Gmats{2};
oT0_Gmat = Gmats{3};
oT1_Gmat = Gmats{4};

end

