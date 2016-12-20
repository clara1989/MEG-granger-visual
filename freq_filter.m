
%% Split into 0.365s non-overlapping epochs
cfg = [];
cfg.channel = data_filtered.label;
cfg.latency = [0.0 0.365];
evoked = ft_selectdata(cfg,data_filtered);
cfg.channel = data_filtered.label;
cfg.latency = [0.365 0.73];
early = ft_selectdata(cfg,data_filtered);
cfg.latency = [0.73 1.095];
middle = ft_selectdata(cfg,data_filtered);
cfg.latency = [1.095 1.46];
late = ft_selectdata(cfg,data_filtered);

cfg = [];
cfg.channel = data_filtered.label;
pstgrating = ft_appenddata(cfg,early,middle,late);

cfg            = [];
cfg.foi = [1:1:140];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.taper      = 'dpss';
cfg.tapsmofrq  = 4;
cfg.keeptrials = 'yes';
cfg.pad = 1;
cfg.padtype = 'zero';
freq    = ft_freqanalysis(cfg, pstgrating);

%% Freq * Spatial Filter

nchan = numel(freq.label);
nfreq  = numel(freq.freq);
nrpttap = size(freq.fourierspctrm,1);

fourier_source = zeros(2,nfreq,nrpttap);

    xxx = spatialfilter(1,:)*reshape(permute(freq.fourierspctrm,[2 3 1]), [nchan nfreq*nrpttap]);
    fourier_source(1,:,:) = reshape(xxx, [nfreq nrpttap]);

freq.fourierspctrm = fourier_source;
freq.dimord = 'chan_freq_rpttap';
freq.label = {'1','2'};

cfg           = [];
cfg.method    = 'granger';
cfg.granger.sfmethod = 'bivariate';
granger      = ft_connectivityanalysis(cfg, freq);

freqflip = ft_freqanalysis(cfg, tmpdata);