%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calculates the peak gamma-band power in the post grating
% period from V1 Virtual Electrode defined using the HCP MMP1.0.. and then
% uses the peak value (40-80Hz) to align granger spectra outputs.
%
% N.B. Make sure to run VE_for_PAC.m or use updated granger_visual_ASD.m
% script to obtain virtsensV1
%
% Written by Robert Seymour January 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%subject = sort({'RS','DB','MP','GW','GR','SY','DS','EC','VS','LA','AE','SW','DK','LH'});
subject = sort({'0401','0402','0403','0404','0405','0406','0407','0409','0411','0413','0414'});  

granger_peak = [];

for i = 1:length(subject)
    
    cd(sprintf('D:\\ASD_Data\\%s\\visual\\PAC',subject{i}))
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
        
        cd(sprintf('D:\\ASD_Data\\%s\\visual\\granger',subject{i}));
        
        % Plot
        collapsed_time = mean(multitaper_diff,3);
        collapsed_time = reshape(collapsed_time,[1,91]);
        subplot(4,4,i);
        x = [30:1:120];plot(x,collapsed_time,'LineWidth',3);
        xlabel('Freq (Hz)');ylabel('Power'); title(subject{i});
        saveas(gcf,'gamma_peak.png')
        
        % Display Maximum Gamma Freq
        max_GBP = 30 + find(collapsed_time == max(collapsed_time(10:end)));
        disp(['Gamma Peak is at ' num2str(max_GBP) 'Hz for subject ' num2str(subject{i}) ]);
        granger_peak(i) = max_GBP;
end

%% Now onto Granger Causality

clear i
%subject = sort({'RS','DB','MP','GW','GR','SY','DS','EC','VS','LA','AE','SW','DK','LH'});
subject = sort({'0401','0402','0403','0405','0406','0409','0411','0413','0414'});  

feedforward = [];
feedback = [];

for i=1:length(subject)
    cd(sprintf('D:\\ASD_Data\\%s\\visual\\granger',subject{i}'));
    load('granger_L.mat'); 
    load('granger_R.mat');
    granger = (granger_L.grangerspctrm + granger_R.grangerspctrm)./2;
    
    ff = [1,3,5,7,9,11,13,15,17,19,21,23,26,27,30];
    fb = [2,4,6,8,10,12,14,16,18,20,22,24,25,28,29];
    
    for j = 1:length(ff)
        feedforward = vertcat(feedforward,granger(ff(j),granger_peak(i)-20:granger_peak(i)+20));
        
    end
for k = 1:length(ff)
        feedback = vertcat(feedback,granger(fb(k),granger_peak(i)-20:granger_peak(i)+20));
    end
end

mean_feedforward = mean(feedforward);
mean_feedback = mean(feedback);
x = [-20:1:20]; figure; plot(x,mean_feedforward,x,mean_feedback)
