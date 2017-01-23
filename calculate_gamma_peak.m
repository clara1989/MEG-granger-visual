%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calculates the peak gamma-band power in the post grating
% period from V1 Virtual Electrode defined using the HCP MMP1.0
%
% N.B. Make sure to run VE_for_PAC.m or use updated granger_visual_ASD.m
% script to obtain virtsensV1
%
% Written by Robert Seymour January 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subject = {'MP','GW','GR','SY','DS','EC','VS','LA','AE'};
subject = {'SW'};

for i = 1:length(subject)
    
    cd(sprintf('D:\\pilot\\%s\\visual\\PAC',subject{i}))
    try
        load('virtsensV1')
    catch
        disp(['Could not load virtsensV1 for subject ' num2str(subject{i})])
    end    
        % TF Analysis
        cfg = [];
        cfg.method = 'mtmconvol';
        cfg.output = 'pow';
        cfg.foi = 30:1:120;
        cfg.toi = [0.3:0.02:1.5];
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;
        cfg.tapsmofrq  = ones(length(cfg.foi),1).*8;
        multitaper_post = ft_freqanalysis(cfg, virtsensV1);
        cfg.toi = [-1.5:0.02:-0.3];
        multitaper_pre = ft_freqanalysis(cfg, virtsensV1);
        
        % Calculate Difference
        multitaper_diff = multitaper_post.powspctrm - multitaper_pre.powspctrm;
        
        cd(sprintf('D:\\pilot\\%s\\visual\\granger',subject{i}));
        
        % Plot
        collapsed_time = mean(multitaper_diff,3);
        collapsed_time = reshape(collapsed_time,[1,91]);
        figure;
        x = [30:1:120];plot(x,collapsed_time,'LineWidth',3);
        xlabel('Freq (Hz)');ylabel('Power'); title(subject{i});
        saveas(gcf,'gamma_peak.png')
        
        % Display Maximum Gamma Freq
        max_GBP = 30 + find(collapsed_time == max(collapsed_time(10:end)));
        disp(['Gamma Peak is at ' num2str(max_GBP) 'Hz for subject ' num2str(subject{i}) ]);
end