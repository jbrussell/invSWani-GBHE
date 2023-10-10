function [ stats, ensemble ] = get_stats( ensemble, par )
% Calculate confidence intervals for results
%
flds = fields(ensemble);
for ifld = 1:length(flds)
    fld = flds{ifld};
    if isstruct(ensemble.(fld))
        jflds = fields(ensemble.(fld));
        for jfld = 1:length(jflds)
%             stats.jflds{jfld}.median = 
            val = ensemble.(fld).(jflds{jfld});
            if strcmp(fld,'E') || strcmp(jflds{jfld}(end),'4')
                I_good = ensemble.X2red_4 < par.chi2red_thresh;
                ensemble.I_good_X2red4 = I_good;
            else
                I_good = ensemble.X2red_2 < par.chi2red_thresh;
                ensemble.I_good_X2red2 = I_good;
            end
            if contains(lower(jflds{jfld}),'dir') || contains(lower(jflds{jfld}),'phi')
                % Circular data (angles)
                [stats.(fld).([jflds{jfld},'_med']), ...
                 stats.(fld).([jflds{jfld},'_u95']), ...
                 stats.(fld).([jflds{jfld},'_l95']), ...
                 stats.(fld).([jflds{jfld},'_u68']), ...
                 stats.(fld).([jflds{jfld},'_l68']) ] = ...
                                get_95_68_prctile_circ(val(:,I_good)',jflds{jfld});
            else 
                [stats.(fld).([jflds{jfld},'_med']), ...
                 stats.(fld).([jflds{jfld},'_u95']), ...
                 stats.(fld).([jflds{jfld},'_l95']), ...
                 stats.(fld).([jflds{jfld},'_u68']), ...
                 stats.(fld).([jflds{jfld},'_l68']) ] = ...
                                get_95_68_prctile(val(:,I_good)');
            end
                        
             [stats.(fld).([jflds{jfld},'_cnts']), ...
              stats.(fld).([jflds{jfld},'_bins'])] = ...
                            hist(val',par.nbins);
        end
    end
end

end


