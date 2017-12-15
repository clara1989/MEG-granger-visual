%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calculates Granger Causality (GC) between 5 visual ROIs 
% (V1,V4,MT,PIT,V7). 

% Data is from my alien task. The ROIs are defined using the HCP MMP atlas 
% and a 3D cortical mesh defined in fs_LR space.
%
% Analysis Steps:
%
% - Perform source analysis across subjsetc specigfic 3D cortical mesh 
% defined in fs_LR space
% - Create parcel specific spatial filter for the 5 ROIs in the left and
% right hemisphere
% - Fourier representation
% - Left multiply the fourier output by the spatial filter for each region
% - Calculate non-parametric Granger Causality using ft_connectivityanalysis
% - Save the left hemisphere and right hemishpere granger spectra
%
% Written by Robert Seymour February 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subject List
subject = {''} 

%% Start Loop
for i=1:length(subject)
    %% Preload the HCP atlas
    % Here we are using the 4k HCP atlas mesh to define visual ROIs from the subject-specific 4k cortical mesh
    load('D:\scripts\MEG-granger-visual\atlas_MSMAll_4k.mat');
    
    %% Load variables required for source analysis
    load(sprintf('D:\\ASD_Data\\%s\\visual\\data_clean_noICA.mat',subject{i})); % non-ICA'd data
    load(sprintf('D:\\ASD_Data\\%s\\visual\\sourceloc\\sens.mat',subject{i}));
    load(sprintf('D:\\ASD_Data\\%s\\visual\\sourceloc\\seg.mat',subject{i}));
    data_filtered = data_clean_noICA; % for book-keeping
    
    %% Set the current directory
    cd(sprintf('D:\\ASD_Data\\%s\\visual\\granger',subject{i}))
    
    %% Set bad channel list - can change to specific channels if necessary
    chans_included = {'MEG', '-MEG0322', '-MEG2542','-MEG0111','-MEG0532'};
    
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
    
    %% Time-Lock Analysis
    cfg = [];
    cfg.channel = chans_included;
    cfg.covariance = 'yes'; % compute the covariance for single trials, then averFEFe
    cfg.preproc.baselinewindow = [-inf 0];  % reapply the baseline correction
    cfg.keeptrials = 'no';
    timelock1 = ft_timelockanalysis(cfg, data_fix);
    
    %% Setup pre-requisites for source localisation
    % Create headmodel
    cfg        = [];
    cfg.method = 'singleshell';
    headmodel  = ft_prepare_headmodel(cfg, seg);
    
    % Load headshape
    headshape = ft_read_headshape(sprintf('D:\\ASD_Data\\raw_alien_data\\rs_asd_%s_aliens_quat_tsss.fif',lower(subject{i})));
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
%     figure; hold on;
%     ft_plot_headshape(headshape)
%     ft_plot_mesh(sourcespace,'facecolor','w','edgecolor',[0.5, 0.5, 0.5],'facealpha',0.1);
%     dataV1 = ft_plot_mesh(sourcemodel_virt.pos(1:8004,:),'vertexcolor','k');
%     ft_plot_sens(sens, 'style', 'black*')
%     set(gcf,'color','w'); drawnow;
    
    %% Perform source analysis across the mesh
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
    
    % I cannot share this at present but the function simply uses SVD to 
    % extract a single spatial filter from each region of the atlas. You
    % can do this using your own code or extract a single spatial filter
    % from each region using the maximum voxel or middle of the parcel etc.
    
    addpath('D:\scripts\MEG-granger-visual');
    [s,p] = mous_lcmv_parcellate(source,timelock1,'parcellation',atlas,'method','parcellation');
    
    %% Split sensor-level into 0.365s non-overlapping epochs (avoiding ERP)
    cfg = [];
    cfg.channel = data_fix.label;
    cfg.latency = [0.365 0.73];
    early = ft_selectdata(cfg,data_fix);
    cfg.latency = [0.73 1.095];
    middle = ft_selectdata(cfg,data_fix);
    cfg.latency = [1.095 1.46];
    late = ft_selectdata(cfg,data_fix);
    
    cfg = [];
    cfg.channel = data_fix.label;
    pstgrating = ft_appenddata(cfg,early,middle,late);
    
    %% Calculate fourier output from 1-140Hz using a Hanning taper
    cfg             = [];
    cfg.foi         = [1:1:140];
    cfg.toi         = 'all';
    cfg.output      = 'fourier';
    cfg.method      = 'mtmfft';
    cfg.taper       = 'hanning';
    cfg.tapsmofrq   = 4; %4Hz smoothing - OK?
    cfg.pad         = 1;
    cfg.keeptrials  = 'yes';
    cfg.padtype     = 'zero';
    freq            = ft_freqanalysis(cfg, pstgrating);
    
    hemisphere = {'L','R'}; % Compute for L then R hemisphere
    
    for hemi = 1:length(hemisphere)
        
        %% Extract VE for LH
        if hemisphere{hemi} == 'L'
            
            % Extract V1 timecourse for later
            [VE_V1] = get_VE_from_parcellated_filter(data_fix,'L_V1_ROI',p);
            % Save virtsensV1 for later PAC calculation
            cd(sprintf('D:\\ASD_Data\\%s\\visual\\PAC',subject{i}))
            save VE_V1 VE_V1
            cd(sprintf('D:\\ASD_Data\\%s\\visual\\granger',subject{i}'));
            clear VE_V1
            
            % Make sure to put extract_fourier_parcellation in your search
            % path
            
            [VE_V1] = extract_fourier_parcellation(freq,'L_V1_ROI',p);
            [VE_V4] = extract_fourier_parcellation(freq,'L_V4_ROI',p);
            [VE_V7] = extract_fourier_parcellation(freq,'L_V7_ROI',p);
            [VE_MT] = extract_fourier_parcellation(freq,'L_MT_ROI',p);
            [VE_PIT] = extract_fourier_parcellation(freq,'L_PIT_ROI',p);
            
        elseif hemisphere{hemi} == 'R'
            
            %% Extract VE for RH
            
            [VE_V1] = extract_fourier_parcellation(freq,'R_V1_ROI',p);
            [VE_V4] = extract_fourier_parcellation(freq,'R_V4_ROI',p);
            [VE_V7] = extract_fourier_parcellation(freq,'R_V7_ROI',p);
            [VE_MT] = extract_fourier_parcellation(freq,'R_MT_ROI',p);
            [VE_PIT] = extract_fourier_parcellation(freq,'R_PIT_ROI',p);
            
        end
        
        %% Create FT structure to hold the fourier outputs
        
        VE_visual = [];
        VE_visual.label = {'V1','V4','V7','MT','PIT'};
        VE_visual.dimord = 'rpttap_chan_freq';
        VE_visual.freq = freq.freq;
        VE_visual.fourierspctrm(:,1,:) = VE_V1;
        VE_visual.fourierspctrm(:,2,:) = VE_V4;
        VE_visual.fourierspctrm(:,3,:) = VE_V7;
        VE_visual.fourierspctrm(:,4,:) = VE_MT;
        VE_visual.fourierspctrm(:,5,:) = VE_PIT;
        VE_visual.cumsumcnt = freq.cumsumcnt;
        VE_visual.cumtapcnt = freq.cumtapcnt;
        VE_visual.trialinfo = freq.trialinfo;
        VE_visual.cfg = freq.cfg;
        
        %% Compute GC + Plot
        cfg           = [];
        cfg.method    = 'granger';
        cfg.granger.sfmethod = 'bivariate';
        granger      = ft_connectivityanalysis(cfg, VE_visual);
        
        cfg = [];
        cfg.parameter = 'grangerspctrm';
        cfg.xlim      = [0 80];
        cfg.zlim = [0 0.05]
        h = figure; ft_connectivityplot(cfg,granger);
        
        if hemisphere{hemi} == 'L'           % if doing LH save accordingly
            saveas(gcf,'granger_L.png')
            granger_L = granger;
            save granger_L granger_L
            VE_visual_L = VE_visual;
            save VE_visual_L VE_visual_L
        end
        
        if hemisphere{hemi} == 'R'           % if doing RH save accordingly
            saveas(gcf,'granger_R.png')
            granger_R = granger;
            save granger_R granger_R
            VE_visual_R = VE_visual;
            save VE_visual_R VE_visual_R
        end
        
    end
    
    %% Display ff versus fb
    
    load('granger_L.mat'); load('granger_R.mat');
    granger = (granger_L.grangerspctrm + granger_R.grangerspctrm)./2;
    granger_L.grangerspctrm = granger;
    
    ff = [1,3,5,7,9,11,13,15,18,19];
    fb = [2,4,6,8,10,12,14,16,17,20];
    
    feedforward = []; feedback = [];
    
    for j = 1:length(ff)
        feedforward = vertcat(feedforward,granger_L.grangerspctrm(ff(j),:));
    end
    for k = 1:length(ff)
        feedback = vertcat(feedback,granger_L.grangerspctrm(fb(k),:));
    end
    
    %% Create a figure
    figure; x = granger_R.freq(1:140);
    plot(x,mean(feedback));hold on; plot(x,mean(feedforward))
    ylim([0 0.06]);xlabel('Frequency (Hz)');ylabel('Granger Causality')
    title([subject{i}]);legend('feedforward','feedback')
    saveas(gcf,'granger_collapsed.png')
end






