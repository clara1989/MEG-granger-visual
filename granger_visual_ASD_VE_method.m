%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calculates Granger Causality (GC) between 6 visual ROIs. Data
% is from my alien task. The ROIs are defined using the HCP MMP atlas and a
% 3D cortical mesh defined in fs_LR space.
%
% Written by Robert Seymour January 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subject List
subject = sort({'RS','DB','MP','GR','DS','EC','VS','LA','AE','SY','GW',...
    'SW','DK','LH','KM','FL','AN','IG'});

%% Start Loop
for i=1:length(subject)
    %% Preload the HCP atlas
    
    % Here we are using the 4k HCP atlas mesh to define visual ROIs from the subject-specific 4k cortical mesh
    load('D:\scripts\MEG-granger-visual\atlas_MSMAll_4k.mat');
    
    %% Load variables required for source analysis
    load(sprintf('D:\\pilot\\%s\\visual\\data_clean_noICA.mat',subject{i}))
    load(sprintf('D:\\pilot\\%s\\visual\\sourceloc\\mri_realigned.mat',subject{i}))
    load(sprintf('D:\\pilot\\%s\\visual\\sourceloc\\sens.mat',subject{i}))
    load(sprintf('D:\\pilot\\%s\\visual\\sourceloc\\seg.mat',subject{i}))
    data_filtered = data_clean_noICA;
    
    %% Set the current directory
    cd(sprintf('D:\\pilot\\%s\\visual\\granger',subject{i}))
    
    %% Set bad channel list - can change to specific channels if necessary
    
    %         addpath('D:\scripts\MEG_preprocessing');
    %         bad_chan_visual_ASD;
    %         chans_included = eval(['chan_list_' subject{i}])
    
    chans_included = {'MEG', '-MEG0322', '-MEG2542'};
    
    %% Load 3D 4k Cortical Mesh for L/R hemisphere & Concatenate
    
    sourcespace = ft_read_headshape({['Subject' subject{i} '.L.midthickness.4k_fs_LR.surf.gii'],['Subject' subject{i} '.R.midthickness.4k_fs_LR.surf.gii']});
    %figure; ft_plot_mesh(sourcespace);camlight; drawnow; %plot if needed
    
    %% Do your timelock analysis on the data & compute covariance
    
    % determine numcomponent by doing an eig on the covariance matrix
    covar = zeros(numel(data_clean_noICA.label));
    for itrial = 1:numel(data_clean_noICA.trial)
        currtrial = data_clean_noICA.trial{itrial};
        covar = covar + currtrial*currtrial.';
    end
    [V, D] = eig(covar);
    D = sort(diag(D),'descend');
    D = D ./ sum(D);
    Dcum = cumsum(D);
    numcomponent = find(Dcum>.99,1,'first'); % number of components accounting for 99% of variance in covar matrix
    
    % Make sure the rank is below 64
    if numcomponent > 65
        numcomponent = 64;
    end
    
    disp(sprintf('\n Reducing the data to %d components \n',numcomponent));
    
    cfg = [];
    cfg.method = 'pca';
    cfg.updatesens = 'yes';
    cfg.channel = chans_included;
    comp = ft_componentanalysis(cfg, data_clean_noICA);
    
    cfg = [];
    cfg.updatesens = 'yes';
    cfg.component = comp.label(numcomponent:end);
    data_fix = ft_rejectcomponent(cfg, comp);
    
    % Time-Lock Analysis
    cfg = [];
    cfg.channel = chans_included;
    cfg.covariance = 'yes'; % compute the covariance for single trials, then averFEFe
    cfg.preproc.baselinewindow = [-inf 0];  % reapply the baseline correction
    cfg.keeptrials = 'no';
    timelock1 = ft_timelockanalysis(cfg, data_fix);
    
    %% Setup pre-requisites for source localisation
    %Create headmodel
    cfg        = [];
    cfg.method = 'singleshell';
    headmodel  = ft_prepare_headmodel(cfg, seg);
    
    % Load headshape
    headshape = ft_read_headshape(sprintf('D:\\pilot\\raw_alien_data\\rs_asd_%s_aliens_quat_tsss.fif',lower(subject{i})));
    headshape = ft_convert_units(headshape,'mm');
    
    %% Create leadfields
    cfg=[];
    cfg.vol=headmodel;
    cfg.channel= chans_included;
    cfg.grid.pos= sourcespace.pos;
    cfg.grid.unit      ='mm';
    cfg.grad=sens;
    cfg.grid.inside = [1:1:length(cfg.grid.pos)]; %always inside - check manually
    cfg.normalize = 'yes';
    sourcemodel_virt=ft_prepare_leadfield(cfg);
    
    % Create Figure to Show Forward Solution
    figure; hold on;
    ft_plot_headshape(headshape)
    ft_plot_mesh(sourcespace,'facecolor','w','edgecolor',[0.5, 0.5, 0.5],'facealpha',0.1);
    dataV1 = ft_plot_mesh(sourcemodel_virt.pos(1:8004,:),'vertexcolor','k');
    ft_plot_sens(sens, 'style', 'black*')
    set(gcf,'color','w'); drawnow;
    
    %% Perform source analysis
    cfg=[];
    cfg.keeptrials = 'no';
    cfg.channel= chans_included;
    cfg.grad = sens;
    cfg.senstype = 'MEG';
    cfg.method='lcmv';
    cfg.grid = sourcemodel_virt;
    cfg.grid.unit      ='mm';
    cfg.headmodel=headmodel;
    cfg.lcmv.lamda='20%';
    cfg.lcmv.fixedori = 'yes';
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.projectmom = 'no';
    cfg.normalize = 'yes';
    source=ft_sourceanalysis(cfg, timelock1);
    
    %% JM Code to get a single spatial filter across all the points of a parcel
    
    [s,p] = mous_lcmv_parcellate(source,timelock1,'parcellation',atlas,'method','parcellation');
    
    %% Now onto the Granger Calculation
    
    hemisphere = {'L','R'};
    
    for hemi = 1:length(hemisphere) %Compute for LH then RH 
        
        %Create figure to show position of VE
        if strcmp(hemisphere{hemi},'L')
            figure; hold on;
            pntV1 = find(atlas.parcellation==2); pntV1 = pntV1(12:49);
            pntV4 = find(atlas.parcellation==7); pntV4 = pntV4(1:13);
            pntMT = find(atlas.parcellation==24);
            pntV7 = find(atlas.parcellation==17);
            pntPIT = find(atlas.parcellation==23);
            
            dataV1 = ft_plot_mesh(sourcemodel_virt.pos(pntV1,:),'vertexcolor','r');
            dataV4 = ft_plot_mesh(sourcemodel_virt.pos(pntV4,:),'vertexcolor','g');
            dataMT = ft_plot_mesh(sourcemodel_virt.pos(pntMT,:),'vertexcolor',[0,0,0]);
            dataV7 = ft_plot_mesh(sourcemodel_virt.pos(pntV7,:),'vertexcolor',[1,1,0]);
            dataPIT = ft_plot_mesh(sourcemodel_virt.pos(pntPIT,:),'vertexcolor',[1,0,1]);
            
            legend([dataV1,dataV4,dataMT,dataV7,dataPIT],'V1','V4','MT','V7','PIT')
            title([subject{i}]);
            hold off
            %legend([dataV1],'V1');% dataV2 dataV4 dataMT dataV7 dataPIT],'V1','V2','V4','MT','V7','PIT');
            ft_plot_mesh(sourcespace,'facecolor','w','edgecolor','none','facealpha',0.1); camlight;
            set(gcf,'color','w'); view(-25, 45); drawnow;
            saveas(gcf,sprintf('VE_loc_%sH.png',(hemisphere{hemi})));
        end
        
        if strcmp(hemisphere{hemi},'R')
            %Create figure to show position of VE
            figure; hold on;
            pntV1 = find(atlas.parcellation==183); pntV1 = pntV1(12:49);
            pntV4 = find(atlas.parcellation==188); pntV4 = pntV4(1:13);
            pntMT = find(atlas.parcellation==205);
            pntV7 = find(atlas.parcellation==198);
            pntPIT = find(atlas.parcellation==204);
            
            dataV1 = ft_plot_mesh(sourcemodel_virt.pos(pntV1,:),'vertexcolor','r');
            dataV4 = ft_plot_mesh(sourcemodel_virt.pos(pntV4,:),'vertexcolor','g');
            dataMT = ft_plot_mesh(sourcemodel_virt.pos(pntMT,:),'vertexcolor',[0,0,0]);
            dataV7 = ft_plot_mesh(sourcemodel_virt.pos(pntV7,:),'vertexcolor',[1,1,0]);
            dataPIT = ft_plot_mesh(sourcemodel_virt.pos(pntPIT,:),'vertexcolor',[1,0,1]);
            
            legend([dataV1,dataV4,dataMT,dataV7,dataPIT],'V1','V4','MT','V7','PIT')
            title([subject{i}]);
            hold off
            %legend([dataV1],'V1');% dataV2 dataV4 dataMT dataV7 dataPIT],'V1','V2','V4','MT','V7','PIT');
            ft_plot_mesh(sourcespace,'facecolor','w','edgecolor','none','facealpha',0.1); camlight;
            set(gcf,'color','w'); view(25, 45); drawnow;
            saveas(gcf,sprintf('VE_loc_%sH.png',(hemisphere{hemi})));
        end
        
        
        %% Extract VE for LH
        if hemisphere{hemi} == 'L'
            
            % Make sure get_VE_from_parcellated_filter is in your search
            % path
            
            [VE_V1] = get_VE_from_parcellated_filter(data_fix,'L_V1_ROI',p);
            [VE_V4] = get_VE_from_parcellated_filter(data_fix,'L_V4_ROI',p);
            [VE_V7] = get_VE_from_parcellated_filter(data_fix,'L_V7_ROI',p);
            [VE_MT] = get_VE_from_parcellated_filter(data_fix,'L_MT_ROI',p);
            [VE_PIT] = get_VE_from_parcellated_filter(data_fix,'L_PIT_ROI',p);
            
            % Save virtsensV1 for later PAC calculation
            cd(sprintf('D:\\pilot\\%s\\visual\\PAC',subject{i}))
            save VE_V1 VE_V1
            cd(sprintf('D:\\pilot\\%s\\visual\\granger',subject{i}'));
            
        elseif hemisphere{hemi} == 'R'
            
            %% Extract VE for RH
            
            [VE_V1] = get_VE_from_parcellated_filter(data_fix,'R_V1_ROI',p);
            [VE_V4] = get_VE_from_parcellated_filter(data_fix,'R_V4_ROI',p);
            [VE_V7] = get_VE_from_parcellated_filter(data_fix,'R_V7_ROI',p);
            [VE_MT] = get_VE_from_parcellated_filter(data_fix,'R_MT_ROI',p);
            [VE_PIT] = get_VE_from_parcellated_filter(data_fix,'R_PIT_ROI',p);
            
        end
        
        %% Append the data & show TFR plots & averaged timeseries
        virtsensparcel=ft_appenddata([],VE_V1,VE_V4,VE_V7,VE_MT,VE_PIT);
        
        cfg = [];
        cfg.channel = virtsensparcel.label;
        cfg.method = 'mtmconvol';
        cfg.output = 'pow';
        cfg.foi = 30:2:120;
        cfg.toi = -2.0:0.02:2.0;
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;
        cfg.tapsmofrq  = ones(length(cfg.foi),1).*8;
        multitaper = ft_freqanalysis(cfg, virtsensparcel);
        
        figure
        for yuy=1:length(multitaper.label)
            cfg=[];
            cfg.channel = multitaper.label{yuy};
            cfg.baselinetype    = 'absolute';
            cfg.ylim            = [40 120];
            cfg.baseline        = [-1.5 0];
            cfg.xlim            = [-1.5 1.5];
            
            subplot(2,3,yuy);ft_singleplotTFR(cfg,multitaper);
            xlabel('Time (sec)'); ylabel('Freq (Hz)');
            colormap(jet)
        end;
        set(gcf, 'Position', get(0,'Screensize')); drawnow;
        
        % Split into 0.365s non-overlapping epochs
        cfg = [];
        cfg.channel = virtsensparcel.label;
        cfg.latency = [0.365 0.73];
        early = ft_selectdata(cfg,virtsensparcel);
        cfg.latency = [0.73 1.095];
        middle = ft_selectdata(cfg,virtsensparcel);
        cfg.latency = [1.095 1.46];
        late = ft_selectdata(cfg,virtsensparcel);
        
        cfg = [];
        cfg.channel = virtsensparcel.label;
        pstgrating = ft_appenddata(cfg,early,middle,late);
        
        %% Calculate granger-causality between all visual areas using a
        %  non-parametric approach
        %  See http://dx.doi.org/10.1016/j.neuron.2015.12.018
        
        % Random time-series for comparison
        tmpdata=virtsensparcel;
        for jj=1:length(tmpdata.trial)
            randtimeseries = randn(size(tmpdata.trial{1,1}));
            sRate = 1000;
            tmpdata.trial{jj}= bst_bandpass_fft(randtimeseries, sRate, 0, sRate/3, 1, 0);
        end;
        clear jj
        
        cfg             = [];
        cfg.foi         = [1:1:140];
        cfg.toi         = pstgrating.time{1, 1};
        cfg.output      = 'fourier';
        cfg.method      = 'mtmfft';
        cfg.taper       = 'hanning';
        cfg.tapsmofrq   = 4;
        cfg.pad         = 4;
        cfg.keeptrials  = 'yes';
        cfg.padtype     = 'zero';
        freq            = ft_freqanalysis(cfg, pstgrating);
        freqflip        = ft_freqanalysis(cfg, tmpdata);
        
        cfg           = [];
        cfg.method    = 'granger';
        cfg.granger.sfmethod = 'bivariate';
        granger      = ft_connectivityanalysis(cfg, freq);
        grangerflip      = ft_connectivityanalysis(cfg, freqflip);
        
        % Show granger spectrum between all 4 areas & save a plot
        cfg = [];
        cfg.parameter = 'grangerspctrm';
        cfg.xlim      = [0 80];
        cfg.zlim = [0 0.05]
        h = figure; ft_connectivityplot(cfg,granger,grangerflip);
        set(h,'NumberTitle','off')           % removes the 'Figure X' title where X is an integer
        set(h,'Name',['GC Spectrum for ', subject{i},' ',hemisphere{hemi},' Hemisphere'])             % set name of window to participant variable
        drawnow;
        
        if hemisphere{hemi} == 'L'           % if doing LH save accordingly
            saveas(gcf,'granger_L.png')
            granger_L = granger;
            save granger_L granger_L
            virtsensparcel_L = virtsensparcel;
            save virtsensparcel_L virtsensparcel_L
        end
        
        if hemisphere{hemi} == 'R'           % if doing RH save accordingly
            saveas(gcf,'granger_R.png')
            granger_R = granger;
            save granger_R granger_R
            virtsensparcel_R = virtsensparcel;
            save virtsensparcel_R virtsensparcel_R
        end
        
    end
    
    %% Using Output Show Peak Granger Spectra
    load('granger_L.mat'); load('granger_R.mat');
    granger = (granger_L.grangerspctrm + granger_R.grangerspctrm)./2;
    granger_L.grangerspctrm = granger;
    
    M = mean(granger)
    r = mean(grangerflip.grangerspctrm)
    
    x = granger_L.freq(1:140); figure; plot(x,M(1:140))
    hold on; plot(x,r(1:140),':');
    xlabel('Frequency (Hz)'); ylabel('Granger Causality')
    legend('Visual','Visual Flipped'); title([subject{i}])
    saveas(gcf,'granger_collapsed.png')
end

