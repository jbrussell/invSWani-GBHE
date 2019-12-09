function [ ensemble ] = collect_models( models, data, par, flds )
% Transform bootstrap output into useful matrices.
%
global data_fds
I.R0 = data.rayl.mode_br_ani == 0;
I.R1 = data.rayl.mode_br_ani == 1;
I.R = logical(I.R0 + I.R1);
I.L0 = data.love.mode_br_ani == 0;
I.L1 = data.love.mode_br_ani == 1;
I.L = logical(I.L0 + I.L1);
I.R0_fds = data_fds.rayl.mode_br_ani == 0;
I.R1_fds = data_fds.rayl.mode_br_ani == 1;
I.R_fds = logical(I.R0_fds + I.R1_fds);
I.L0_fds = data_fds.love.mode_br_ani == 0;
I.L1_fds = data_fds.love.mode_br_ani == 1;
I.L_fds = logical(I.L0_fds + I.L1_fds);
for ifld = 1:length(flds)
    ensemble.(flds{ifld}).([flds{ifld},'_mag']) = zeros(size(models(1).(flds{ifld}).z,1),size(models,2));
    ensemble.(flds{ifld}).([flds{ifld},'_dir']) = zeros(size(models(1).(flds{ifld}).z,1),size(models,2));
end
ensemble.X2red_2 = zeros(size(models,2),1);
ensemble.X2red_4 = zeros(size(models,2),1);
iflds = {'R0','R1','R','L0','L1','L','R0_fds','R1_fds','R_fds','L0_fds','L1_fds','L_fds'};
type =  {'rayl','rayl','rayl','love','love','love','rayl','rayl','rayl','love','love','love'};
for ifld = 1:length(iflds)
    jflds = {'A2','phi2','A4','phi4'};        
    for jfld = 1:length(jflds)
        if isfield(models(1).(iflds{ifld}),jflds{jfld})
            ensemble.(iflds{ifld}).(jflds{jfld}) = nan(size(models(1).(iflds{ifld}).(jflds{jfld}),1),size(models,2));
        end
    end
end
    
for ibs = 1:length(models)
    model = models(ibs);
    for ifld = 1:length(flds)
        fld = flds{ifld};
        ensemble.(fld).([fld,'_mag'])(:,ibs) = model.(fld).([fld,'_mag']);
        ensemble.(fld).([fld,'_dir'])(:,ibs) = model.(fld).([fld,'_dir']);
    end
    ensemble.X2red_2(ibs) = model.chi2.X2red_2;
    ensemble.X2red_4(ibs) = model.chi2.X2red_4;
    for ifld = 1:length(iflds)
        for jfld = 1:length(jflds)
            if isfield(model.(iflds{ifld}),jflds{jfld})
                model_pers = model.(iflds{ifld}).periods;
                data_pers = data.(type{ifld}).periods_ani(I.(iflds{ifld}));
                [~,I_pers_m] = intersect(model_pers,data_pers);
                [~,I_pers_d] = intersect(data_pers,model_pers);
                ensemble.(iflds{ifld}).(jflds{jfld})(I_pers_d,ibs) = model.(iflds{ifld}).(jflds{jfld})(I_pers_m);
            end
        end
    end
end
% ensemble.I_good_X2red2(1,:) = ensemble.X2red_2 < par.chi2red_thresh;
% ensemble.I_good_X2red4(1,:) = ensemble.X2red_4 < par.chi2red_thresh;


end

