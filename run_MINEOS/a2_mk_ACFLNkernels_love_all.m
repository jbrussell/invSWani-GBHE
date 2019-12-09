%% Drive to calculate azimuthal anisotropy kernels (same as A,C,F,L,N kernels 
% from Dziewonski & Anderson 1981)
%
% RAYLEIGH WAVES (sensitive to G, B, H):
% "G" kernels for 2-theta SV anisotropy: use "L" kernels
%       "B" kernels for 2-theta SV and PH anisotropy: use "A" kernels (not well constrained by Rayleigh waves)
%       "H" kernels for 2-theta SV: use "F" kernels (not well constraine by Rayleigh waves)
%
% LOVE WAVES (sensitive to G, E):
% "G" kernels for 2-theta SH anisotropy: use "L" kernels
% "E" kernels for 4-theta SH anisotropy: use "N" kernels
%
%
% NOTE: 
% These kernels contain no scaling and therefore describe perturbations to 
% eigenfrequency rather than phase velocity. To get in units of phase velocity
% multiply by c/U. Then, to invert percent perturbations to phase velocity for
% percent perturbations to SV and SH (i.e. G/L, B/A, H/F, E/N), need to multiply
% by rho * V^2. Therefore in order to use these kernels for inversion, need to 
% multiply them by:
%                   c/U * rho * V^2
%
% FIRST RUN MINEOS TO GENERATE MODE TABLES !
% (Should have all files for all desired mode and branches in the ./MODE/TABLES/CARDID)
%
clear
%plot native
SONLY_vec = [1 1 0 0]; %1 spheroidal (0 toroidal)
TONLY_vec = [0 0 1 1]; %1 Toroidal (0 spheroidal)
branch_vec = [0 1 0 1]; %0 fundamental, 1 first overtone

% SONLY_vec = [1 1 0]; %1 spheroidal (0 toroidal)
% TONLY_vec = [0 0 1]; %1 Toroidal (0 spheroidal)
% branch_vec = [0 1 0]; %0 fundamental, 1 first overtone
for imode = 1:length(SONLY_vec)
    clear h lgd Gmatrix
    SONLY = SONLY_vec(imode);
    TONLY = TONLY_vec(imode);
    branch = branch_vec(imode);
    
    % write file for parameter_ACFLN.m to read
    fid = fopen('mode.txt','w');
    fprintf(fid,'%.0f %.0f %.0f',SONLY,TONLY,branch);
    fclose(fid);

    parameter_FRECHET;

    model_depth = 400; %80; %80 km
    issave_Gmatrix = 1; % save G matrix?
    issave_FRECH = 1; % save frechet kernels?

    %TYPE = 'S'
    TYPE = param.TYPE
    % branch = 0;

    CARDID = param.CARDID;
    periods = param.periods;
    FRECHETPATH = param.frechetpath;
    ylims = [0 100]; %[0 400];
    %ylims = [0 700];
    xlims = [0 1e-17]*1000*4;



    isfigure = 1;

    %% Set path to executables
    setpath_plotwk;

    %% Change environment variables to deal with gfortran
    setenv('GFORTRAN_STDIN_UNIT', '5') 
    setenv('GFORTRAN_STDOUT_UNIT', '6') 
    setenv('GFORTRAN_STDERR_UNIT', '0')


    %% run plot_wk on the table_hdr file to generate the branch file
    write_plotwk(TYPE,CARDID);

    com = ['cat run_plotwk.',lower(TYPE),' | plot_wk'];
    [status,log] = system(com);

    if status ~= 0     
        error( 'something is wrong at plot_wk')
    end

    %% run "frechet" to generate the frechet file 
    NDISC = 0;
    ZDISC = [];
    if ( TYPE == 'T') 
        TYPEID = param.TTYPEID;
    elseif ( TYPE == 'S') 
        TYPEID = param.STYPEID;
    end

    com = ['ls ',param.TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_1.eig_fix | cat'];
    [status eig_fils] = system(com);
    if strcmp(eig_fils(end-25:end-1),'No such file or directory')
        disp('Found no *.eig_fix files')
        write_frechet(TYPE,CARDID,NDISC,ZDISC)
    else
        disp('Found *.eig_fix files')
        write_frech_chk(NDISC)
    end
    disp('Be patient! This will take ~25 s');
    tic
    com = ['cat run_frechet.',lower(TYPE),' | frechet_ACFLN_love > frechet.LOG'];
    [status,log] = system(com);
    if status ~= 0     
        error( 'something is wrong at frechet')
    end
    toc

    %% run "frechet_cv" to generate the cv kernels
    % Convert frechet to ascii
    % Make CV Frechet Kernels
    disp('--- Make CV Frechet Kernels ---');


    write_frechcv(TYPE,CARDID,branch)

    com = ['cat run_frechcv.',lower(TYPE),' | frechet_cvG'];
    [status,log] = system(com);
    if status ~= 0     
        error( 'something is wrong at frechet_cvG')
    end
    %% Convert CV Frechet kernels to ascii with phase-velocity sensitivity
    % Will do this for all periods of interest
    % Set inside the setparam_MINE.m
    disp('--- Convert Frechet CV to ascii ---');

        % Program writes run file for draw_frechet_gv, runs it, and reads in
        % sensitivity kernels for all periods of interest

    FRECH = frechACFLN_asc(TYPE,CARDID,branch);  

    %% Save kernels in G matrix
    Idep = find(FRECH(1).rad >= FRECH(1).rad(end)-model_depth*1000);
    for ip = 1:length(periods)
        if ( TYPE == 'T') 
            Gmatrix.L(ip,:) = flipud(FRECH(ip).L(Idep));
            Gmatrix.N(ip,:) = flipud(FRECH(ip).N(Idep)); 
        elseif ( TYPE == 'S') 
            Gmatrix.A(ip,:) = flipud(FRECH(ip).A(Idep));
            Gmatrix.C(ip,:) = flipud(FRECH(ip).C(Idep));
            Gmatrix.F(ip,:) = flipud(FRECH(ip).F(Idep));
            Gmatrix.L(ip,:) = flipud(FRECH(ip).L(Idep));
        end
    end
    Gmatrix.rad = flipud(FRECH(ip).rad(Idep));
    Gmatrix.model_depth = model_depth;
    Gmatrix.periods = periods;
    depthlayer = (6371000-Gmatrix.rad)/1000;
    if issave_Gmatrix
    %     save([FRECHETPATH,'Gmatrix_',TYPE,num2str(branch),'_',CARDID,'.mat'],'Gmatrix','depthlayer');
        save([FRECHETPATH,'Gmatrix_',TYPE,num2str(branch),'_',num2str(periods(1)),'_',num2str(periods(end)),'s_',CARDID,'.mat'],'Gmatrix','depthlayer');
    end

    %% Plot Kernels

    fig1 = figure(1); clf;
    %set(gcf,'position',[58   255   548   450]);
    % set(gcf,'position',[58   255   916   450]);

    clr = jet(length(periods));
    rad = (FRECH(1).rad(end)-FRECH(1).rad)/1000;               
    for iper = 1:length(FRECH) %length(periods)
        lgd{iper} = [num2str(FRECH(iper).per),' s'];

        if ( TYPE == 'T') 
            set(gcf,'position',[58   255   916   450],'color','w');

            % L
            subplot(1,3,1);
            set(gca,'linewidth',2);
            h(iper) = plot(FRECH(iper).L,rad,'linewidth',2,'color',clr(iper,:)); hold on;
            axis ij
            title(['L (2\theta)'],'fontsize',15);
            ylabel('Depth (km)','fontsize',15);
            ylim(ylims);
            %xlim([0 2e-5]);
            xlim(xlims*30);
            set(gca,'fontsize',15)

            % N
            subplot(1,3,2);
            set(gca,'linewidth',2);
            plot(FRECH(iper).N,rad,'linewidth',2,'color',clr(iper,:)); hold on;
            axis ij
            title(['N (4\theta)'],'fontsize',15);
            ylabel('Depth (km)','fontsize',15);
            ylim(ylims);
    %         xlim([1e-16 1e-11]);
            xlim(xlims*200);
            set(gca,'fontsize',15)

            % RHO
            subplot(1,3,3);
            set(gca,'linewidth',2);
            plot(FRECH(iper).rho,rad,'linewidth',2,'color',clr(iper,:)); hold on;
            axis ij
            title(['\rho Kernels'],'fontsize',15);
            ylabel('Depth (km)','fontsize',15);
            ylim(ylims);
    %         xlim([1e-16 1e-11]);
            xlim(xlims);
            set(gca,'fontsize',15)



        elseif ( TYPE == 'S') 
            set(gcf,'position',[10         234        1174         471]);

            % A
            subplot(1,5,1);
            set(gca,'linewidth',2);
            h(iper) = plot(FRECH(iper).A,rad,'linewidth',2,'color',clr(iper,:)); hold on;
            axis ij
            title(['A (2\theta & 4\theta)'],'fontsize',15);
            ylabel('Depth (km)','fontsize',15);
            ylim(ylims);
    %         xlim([1e-16 1e-11]);
            xlim(xlims*15);
            set(gca,'fontsize',15)

            % C
            subplot(1,5,2);
            set(gca,'linewidth',2);
            plot(FRECH(iper).C,rad,'linewidth',2,'color',clr(iper,:)); hold on;
            axis ij
            title(['C'],'fontsize',15);
            ylabel('Depth (km)','fontsize',15);
            ylim(ylims);
    %         xlim([1e-16 1e-11]);
            xlim(xlims*5);
            set(gca,'fontsize',15)

            % F
            subplot(1,5,3);
            set(gca,'linewidth',2);
            if iper == 1
                plot([0 0],ylims,'--k'); hold on;
            end
            plot(FRECH(iper).F,rad,'linewidth',2,'color',clr(iper,:)); hold on;
            axis ij
            title(['F (2\theta)'],'fontsize',15);
            ylabel('Depth (km)','fontsize',15);
            ylim(ylims);
    %         xlim([1e-16 1e-11]);
            xlim(xlims*5);
            set(gca,'fontsize',15)

            % L
            subplot(1,5,4);
            set(gca,'linewidth',2);
            plot(FRECH(iper).L,rad,'linewidth',2,'color',clr(iper,:)); hold on;
            axis ij
            title(['L (2\theta)'],'fontsize',15);
            ylabel('Depth (km)','fontsize',15);
            ylim(ylims);
    %         xlim([1e-16 1e-11]);
            xlim(xlims*20);
            set(gca,'fontsize',15)

            % RHO
            subplot(1,5,5);
            set(gca,'linewidth',2);
            plot(FRECH(iper).rho,rad,'linewidth',2,'color',clr(iper,:)); hold on;
            axis ij
            title(['\rho'],'fontsize',15);
            ylabel('Depth (km)','fontsize',15);
            ylim(ylims);
    %         xlim([1e-16 1e-11]);
            xlim(xlims);
            set(gca,'fontsize',15)
        end
        %display(FRECH(iper).per)
        %pause;
    end

    if ( TYPE == 'T')
        subplot(1,3,1); hold on;
        legend(h,lgd,'location','southeast','fontsize',15);
    elseif ( TYPE == 'S') 
        subplot(1,5,1); hold on;
        legend(h,lgd,'location','southeast','fontsize',15);
        % plot([0,0],[0,3000],'--k');
    end


    % subplot(1,2,2); hold on;
    % legend(lgd,'location','southeast');
    % plot([0,0],[0,3000],'--k');
    % get(h,'interpreter');
    % set(h,'interpreter','none');

    CARDID = param.CARDID;
    % TYPEID = param.TTYPEID;
    %print('-painters','-dpdf','-r400',[EIGPATH,CARDID,'.',TYPEID,'.',num2str(j),'mod.',num2str(N_modes),'_fix.pdf']);
    save2pdf([FRECHETPATH,'ACFLNkernels_',TYPE,'_',num2str(periods(1)),'_',num2str(periods(end)),'_',num2str(ylims(2)),'km.',num2str(branch),'.pdf'],fig1,1000);

    if exist(['run_plotwk.',lower(TYPE)])
        delete(['run_plotwk.',lower(TYPE)])
    end
    delete(['run_frechcv.',lower(TYPE)],['run_frechet.',lower(TYPE)],['run_frechcv_asc.',lower(TYPE)]);
    % savefile = [FRECHETPATH,'fACFLN_',TYPE,num2str(branch),'_',CARDID,'.mat'];
    savefile = [FRECHETPATH,'fACFLN_',TYPE,num2str(branch),'_',num2str(periods(1)),'_',num2str(periods(end)),'s_',CARDID,'.mat'];
    if issave_FRECH
        save(savefile,'FRECH');
    end

    % Copy *.q file to FRECH directory
    QPATH = [param.TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.q'];
%     if ~exist([FRECHETPATH,CARDID,'.',TYPEID,'.q'])
    system(['cp ',QPATH,' ',FRECHETPATH ])
%     end
    
    % Delete large frechet files
    delete([FRECHETPATH,CARDID,'.',TYPEID,'.fcv.',num2str(branch)])
    delete([FRECHETPATH,CARDID,'.',TYPEID,'.frech'])
end

% Change the environment variables back to the way they were
setenv('GFORTRAN_STDIN_UNIT', '-1') 
setenv('GFORTRAN_STDOUT_UNIT', '-1') 
setenv('GFORTRAN_STDERR_UNIT', '-1')

delete mode.txt frechet.LOG draw_frechet_gv.LOG
