%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Granger Causality Statistical Analysis
%
% Need to figure out a way of exporting the mask for pos + neg clusters
% separately so that it can be plotted alongside my granger spectra
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subject = sort({'RS','DB','MP','GR','DS','EC','VS','LA','AE','SY','GW',...
%     'SW','DK','LH','KM','FL','AN','IG'});

subject = {'0401','0402','0403','0404','0405','0406','0407','0409','0411',...
     '0413','0414','0415','0416','0417'};
 
grandavgA = []

for i=1:length(subject)
    cd(sprintf('D:\\ASD_Data\\%s\\visual\\granger',subject{i}'));
    % Load in Data
    load('granger_L.mat');load('granger_R.mat');
    % Average Over Hemispheres - can also change + run individually
    granger = (granger_L.grangerspctrm + granger_R.grangerspctrm)./2;
    
    % Put into FT Structure with Freq = Time
    granger_ASD_Data = [];
    granger_ASD_Data.label = {'ff'};
    granger_ASD_Data.dimord = 'chan_time';
    granger_ASD_Data.avg(1,:) = granger(1,:);
    
    granger_ASD_Data.time = [1:1:140];
    % Add to grandavgA
    grandavgA{i} = granger_ASD_Data;
end

grandavgB = [];

for i=1:length(subject)
    cd(sprintf('D:\\ASD_Data\\%s\\visual\\granger',subject{i}'));
    % Load in Data
    load('granger_L.mat');load('granger_R.mat');
    % Average Over Hemispheres
    granger = (granger_L.grangerspctrm + granger_R.grangerspctrm)./2;
    
    % Put into FT Structure with Freq = Time
    granger_ASD_Data = [];
    granger_ASD_Data.label = {'ff'};
    granger_ASD_Data.dimord = 'chan_time';
    granger_ASD_Data.avg(1,:) = granger(2,:);
    
    granger_ASD_Data.time = [1:1:140];
    % Add to grandavgB
    grandavgB{i} = granger_ASD_Data;
end

%% Cluster-Based Statistics (ft_timelockstatistics) like you would use for ERPs
cfg = [];
cfg.avgovertime = 'no';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.clusteralpha = 0.05;
cfg.correctm    = 'cluster';
cfg.numrandomization = 1000;
 
Nsub = numel(grandavgA);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,grandavgB{:},grandavgA{:});

figure; plot(stat.stat); xlabel('Freq (Hz)'); ylabel('t-value');
figure; plot(stat.prob);xlabel('Freq (Hz)'); ylabel('p-value');