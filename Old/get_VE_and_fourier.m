%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate virtual electrode time-series and source-level
% fourier representions.
%
% The script uses PCA to define a single spatial filter from a particular
% region and then multiples this by senso-level data (VE) and the
% sensor-level fourier output

% Inputs:
% - source = output of ft_sourceanlysis (LCMV beamformer) with
% cfg.keepfilter = 'yes' and cfg.fixedori = 'yes'
% - indices = indices of the grid points of interest
% - freq = sensor level fourier freq structure
% - avg = sensor-level covrince matrix (ft_timelockanalysis)
% - label_atlas = name of the region
% - data_clean = sensor-level data
%
% Outputs:
% - VE = virtual elctrode time-series
% - fourier = sensor-level fouier-output * spatial filter
%
% Written by Robert Seymour (17/07/17)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [VE,fourier] = get_VE_and_fourier(source,indices,freq,avg,label_atlas, data_clean)

% USe PCA to get spatial filter
F   = cat(1,source.avg.filter{indices(1):indices(2)});
[u,s,v] = svd(F*avg.cov*F');
filter = u'*F;

% Multiply spatial filter by sensor-level fourier output
nchan = numel(freq.label);
nfreq  = numel(freq.freq);
nrpttap = size(freq.fourierspctrm,1);

fourier_source = filter(1,:)*reshape(permute(freq.fourierspctrm,[2 3 1]), [nchan nfreq*nrpttap]);
fourier_source = reshape(fourier_source, [nfreq nrpttap]);
fourier = rot90(fourier_source);

%% This extracts timeseries info for each trial rather than the average
VE = [];
VE.label = {label_atlas};
VE.trialinfo = data_clean.trialinfo;
for d=1:(length(data_clean.trialinfo))
    % note that this is the non-filtered raw data
    VE.time{d}       = data_clean.time{d};
    VE.trial{d}(1,:) = filter(1,:) * data_clean.trial{d}(:,:);
end

%% Compute multitaper for the virtual channel timeseries
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.foi = 30:1:100;
cfg.toi = -2.0:0.02:2.0;
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.25;
cfg.tapsmofrq  = ones(length(cfg.foi),1).*8;
multitaper = ft_freqanalysis(cfg, VE);

%% Plot
cfg                 = [];
%cfg.baselinetype    = 'relative';
cfg.ylim            = [40 100];
cfg.baseline        = [-1.5 0];
cfg.xlim            = [-0.5 1.5]
%cfg.channel         = {'pos'};
figure; ft_singleplotTFR(cfg, multitaper);
colormap(jet);
title(label_atlas)
end
