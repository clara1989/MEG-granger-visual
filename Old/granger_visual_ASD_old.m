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

hemisphere = {'L','R'};

%% Start Loop
for i=1:length(subject)
    %% Preload the HCP atlas    
    
    % Here we are using the 4k HCP atlas mesh to define visual ROIs from the subject-specific 4k cortical mesh
    cifti_atlas = ft_read_cifti('D:\HCP_atlas\HCP_atlas_downsampled_4k.dlabel.nii');
    
    for hemi = 1:2
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
        
        %% Get the Atlas Points of Interest on the HCP atlas
                
        % Get ROIs based on L/R hemisphere
        
        if hemisphere{hemi} == 'L'           % if doing LH these are the atlas labels
            V1 = 1; V2 = 4;V3 = 5; V4 = 6; MT = 23; PIT = 22; V7 = 16;
        end
        
        if hemisphere{hemi} == 'R'           % if doing RH these are the atlas labels
            V1 = 181; V2 = 184;V3 = 185; V4 = 186; MT = 203; PIT = 202; V7 = 196;
        end
        
        % V1
        
        indxV1 = find(cifti_atlas.x1 == V1); indxV1= indxV1(12:49);
        indxV1pnt = zeros(length(indxV1),3);
        clear k;
        for k = 1:length(indxV1)
            indxV1pnt(k,:) = sourcespace.pos(indxV1(k),:);
        end
        
%         % V2
%         
%         indxV2 = find(cifti_atlas.x1 == V2); indxV2 = indxV2(8:14);
%         indxV2pnt = zeros(length(indxV2),3);
%         clear k;
%         for k = 1:length(indxV2)
%             indxV2pnt(k,:) = sourcespace.pos(indxV2(k),:);
%         end
%         
        % V4
        
        indxV4 = find(cifti_atlas.x1 == V4); 
        indxV4pnt = zeros(length(indxV4),3);
        clear k;
        for k = 1:length(indxV4)
            indxV4pnt(k,:) = sourcespace.pos(indxV4(k),:);
        end
        
%         % MT
%         
%         indxMT = find(cifti_atlas.x1 == MT);
%         indxMTpnt = zeros(length(indxMT),3);
%         clear k;
%         for k = 1:length(indxMT)
%             indxMTpnt(k,:) = sourcespace.pos(indxMT(k),:);
%         end
%         
%         % PIT
%         
%         indxPIT = find(cifti_atlas.x1 == PIT);
%         indxPITpnt = zeros(length(indxPIT),3);
%         clear k;
%         for k = 1:length(indxPIT)
%             indxPITpnt(k,:) = sourcespace.pos(indxPIT(k),:);
%         end
%         
%         % V7
%         
%         indxV7 = find(cifti_atlas.x1 == V7);
%         indxV7pnt = zeros(length(indxV7),3);
%         clear k;
%         for k = 1:length(indxV7)
%             indxV7pnt(k,:) = sourcespace.pos(indxV7(k),:);
%         end
               
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
        cfg.keeptrials = 'yes';
        timelock1 = ft_timelockanalysis(cfg, data_fix);
        
        %% Setup pre-requisites for source localisation
        %Create headmodel
        cfg        = [];
        cfg.method = 'singleshell';
        headmodel  = ft_prepare_headmodel(cfg, seg);
        
        % Load headshape
        headshape = ft_read_headshape(sprintf('D:\\pilot\\raw_alien_data\\rs_asd_%s_aliens_quat_tsss.fif',lower(subject{i})));
        headshape = ft_convert_units(headshape,'mm');
        
        %% Create leadfields for visual areas
        cfg=[];
        cfg.vol=headmodel;
        cfg.channel= chans_included;
        cfg.grid.pos=[indxV1pnt; indxV4pnt];% indxMTpnt; indxV7pnt; indxPITpnt];% indxV4pnt; 
        cfg.grid.unit      ='mm';
        cfg.grad=sens;
        cfg.grid.inside = [1:1:length(cfg.grid.pos)]; %always inside - check manually
        cfg.normalize = 'yes';
        sourcemodel_virt=ft_prepare_leadfield(cfg);
        
        % Create Figure to Show Forward Solution
%         figure; hold on;
%         ft_plot_headshape(headshape)
%         ft_plot_mesh(sourcemodel_virt.pos(sourcemodel_virt.inside,:),'white','none');
%         ft_plot_mesh(sourcespace,'facecolor','w','edgecolor',[0.5, 0.5, 0.5],'facealpha',0.1);
%         ft_plot_sens(sens, 'style', 'black*')
%         title([subject{i}]); drawnow;
        
        %% Create grid point diagram for LH
        if strcmp(hemisphere{hemi},'L')
            %Create figure to show position of VE
            figure; hold on;
            dataV1 = ft_plot_mesh(sourcemodel_virt.pos(1:38,:),'vertexcolor','r');
            %       dataV2 = ft_plot_mesh(sourcemodel_virt.pos(12:19,:),'vertexcolor','b');
            dataV4 = ft_plot_mesh(sourcemodel_virt.pos(39:73,:),'vertexcolor','g');
%             dataMT = ft_plot_mesh(sourcemodel_virt.pos(39:45,:),'vertexcolor',[0,0,0]);
%             dataV7 = ft_plot_mesh(sourcemodel_virt.pos(46:51,:),'vertexcolor',[1,1,0]);
%             dataPIT = ft_plot_mesh(sourcemodel_virt.pos(52:60,:),'vertexcolor',[1,0,1]);
            legend([dataV1; dataV4],'V1','V4')
            title([subject{i}]);
            hold off
            %legend([dataV1],'V1');% dataV2 dataV4 dataMT dataV7 dataPIT],'V1','V2','V4','MT','V7','PIT');
            ft_plot_mesh(sourcespace,'facecolor','w','edgecolor','none','facealpha',0.1); camlight;
            set(gcf,'color','w'); view(-25, 45); drawnow;
            saveas(gcf,sprintf('VE_loc_%sH.png',(hemisphere{hemi})));
        end
        
        %% Create grid point diagram for RH
        if strcmp(hemisphere{hemi},'R')
            %Create figure to show position of VE
            figure; hold on;
            dataV1 = ft_plot_mesh(sourcemodel_virt.pos(1:38,:),'vertexcolor','r');
            %       dataV2 = ft_plot_mesh(sourcemodel_virt.pos(12:19,:),'vertexcolor','b');
            dataV4 = ft_plot_mesh(sourcemodel_virt.pos(39:73,:),'vertexcolor','g');
%             dataMT = ft_plot_mesh(sourcemodel_virt.pos(39:46,:),'vertexcolor',[0,0,0]);
%             dataV7 = ft_plot_mesh(sourcemodel_virt.pos(46:52,:),'vertexcolor',[1,1,0]);
%             dataPIT = ft_plot_mesh(sourcemodel_virt.pos(53:58,:),'vertexcolor',[1,0,1]);
            legend([dataV1; dataV4]);%; dataV7; dataPIT],'V1','MT','V7','PIT')
            title([subject{i}]);
            hold off
            %legend([dataV1],'V1');% dataV2 dataV4 dataMT dataV7 dataPIT],'V1','V2','V4','MT','V7','PIT');
            ft_plot_mesh(sourcespace,'facecolor','w','edgecolor','none','facealpha',0.1); camlight;
            set(gcf,'color','w'); view(25, 45); drawnow;
            saveas(gcf,sprintf('VE_loc_%sH.png',(hemisphere{hemi})));
        end
        
        %% Perform source analysis for visual areas
        cfg=[];
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
        
        cfg = [];
        cfg.channel= chans_included;
        data_filtered = ft_preprocessing(cfg,data_fix);
        
        figure
        plot(source.time, source.avg.mom{1})
        
        %% Extract VE timeseries from visual areas
        spatialfilter=cat(1,source.avg.filter{:});
        virtsens=[];
        for nnn=1:length(data_filtered.trial)
            virtsens.trial{nnn}=spatialfilter*data_filtered.trial{nnn};
        end;
        
        virtsens.time=data_filtered.time;
        virtsens.fsample=data_filtered.fsample;
        indx=[indxV1pnt; indxV4pnt];% indxV7pnt;indxPITpnt];
        for nnn=1:length(virtsens.trial{1}(:,1))
            virtsens.label{nnn}=[num2str(nnn)];
        end;
        
        %% Extract each gridpoint and average over the corresponding area
        cfg = [];
        chan_number = length(indxV1);
        cfg.channel = virtsens.label(1:chan_number);
        cfg.avgoverchan = 'yes';
        virtsensV1 = ft_selectdata(cfg,virtsens);
        virtsensV1.label = {'V1'};
        
        % Save virtsensV1 for later PAC calculation
        cd(sprintf('D:\\pilot\\%s\\visual\\PAC',subject{i}))
        save virtsensV1 virtsensV1
        cd(sprintf('D:\\pilot\\%s\\visual\\granger',subject{i}'));
        
%         cfg.channel = virtsens.label(chan_number+1:chan_number+(length(indxV2)))
%         chan_number = chan_number+(length(indxV2));
%         cfg.avgoverchan = 'yes';
%         virtsensV2 = ft_selectdata(cfg,virtsens);
%         virtsensV2.label = {'V2'};
%           
        cfg.channel = virtsens.label(chan_number+1:chan_number+(length(indxV4)))
        chan_number = chan_number+(length(indxV4));
        cfg.avgoverchan = 'yes';
        virtsensV4 = ft_selectdata(cfg,virtsens);
        virtsensV4.label = {'V4'};
        
% %         cfg.channel = virtsens.label(chan_number+1:chan_number+(length(indxMT)))
% %         chan_number = chan_number+(length(indxMT));
% %         cfg.avgoverchan = 'yes';
% %         virtsensMT = ft_selectdata(cfg,virtsens);
% %         virtsensMT.label = {'MT'};
% %         
% %         cfg.channel = virtsens.label(chan_number+1:chan_number+(length(indxV7)))
% %         chan_number = chan_number+(length(indxV7));
% %         cfg.avgoverchan = 'yes';
% %         virtsensV7 = ft_selectdata(cfg,virtsens);
% %         virtsensV7.label = {'V7'};
% %         
% %         cfg.channel = virtsens.label(chan_number+1:chan_number+(length(indxPIT)))
% %         chan_number = chan_number+(length(indxPIT));
% %         cfg.avgoverchan = 'yes';
% %         virtsensPIT = ft_selectdata(cfg,virtsens);
% %         virtsensPIT.label = {'PIT'};
        
        
        %% Append the data & show TFR plots & averaged timeseries
        virtsensparcel=ft_appenddata([],virtsensV1,virtsensV4);%,virtsensV7, virtsensPIT);
        
        cfg = [];
        cfg.channel = virtsensparcel.label;
        cfg.method = 'mtmconvol';
        cfg.output = 'pow';
        cfg.foi = 30:2:120;
        cfg.toi = -2.0:0.02:2.0;
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;
        cfg.tapsmofrq  = ones(length(cfg.foi),1).*8;
        multitaper = ft_freqanalysis(cfg, virtsensparcel);
        
        figure;
        for yuy=1:length(multitaper.label)
            cfg=[];
            cfg.channel = multitaper.label{yuy};
            
            cfg.baselinetype    = 'absolute';
            cfg.ylim            = [40 120];
            cfg.baseline        = [-1.5 0];
            cfg.xlim            = [-1.5 1.5];
            
            subplot(1,2,yuy);ft_singleplotTFR(cfg,multitaper);
            xlabel('Time (sec)'); ylabel('Freq (Hz)');
            colormap(jet)
        end;
        %set(gcf, 'Position', get(0,'Screensize')); drawnow;
                
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
        
        % tmpdata=pstgrating;
        % for jj=1:length(tmpdata.trial)
        %     tmpdata.trial{jj}=fliplr(tmpdata.trial{jj});
        % end;
        
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
        
        % Save the granger spectrum & virtual sensor data
        if hemisphere{hemi} == 'L'
            granger_L = granger; save('granger_L', 'granger_L','-v7.3');
            save virtsensparcel_L virtsensparcel
        end
        
        if hemisphere{hemi} == 'R'
            granger_R = granger; save('granger_R', 'granger_R','-v7.3');
            save virtsensparcel_R virtsensparcel
        end
        
        % Show granger spectrum between all 4 areas & save a plot
        cfg = [];
        cfg.channel = {'V1','V4'};
        cfg.parameter = 'grangerspctrm';
        cfg.xlim      = [0 80];
        cfg.zlim = [0 0.05]
        h = figure; ft_connectivityplot(cfg,granger,grangerflip);
        set(h,'NumberTitle','off')           % removes the 'Figure X' title where X is an integer
        set(h,'Name',['GC Spectrum for ', subject{i},' ',hemisphere{hemi},' Hemisphere'])             % set name of window to participant variable
        drawnow;
        
        if hemisphere{hemi} == 'L'           % if doing LH save accordingly
            saveas(gcf,'granger_L.png')
        end
        
        if hemisphere{hemi} == 'R'           % if doing RH save accordingly
            saveas(gcf,'granger_R.png')   
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

