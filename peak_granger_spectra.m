subject = {'RS','DB','MP','GW','GR','SY','DS','AE','EC'};
granger_comb = [];   
for i=1:length(subject)
    cd(sprintf('D:\\pilot\\%s\\visual\\granger',subject{i}'));
    load('granger_L.mat'); load('granger_R.mat');
    granger = (granger_L.grangerspctrm);% + granger_R.grangerspctrm);
    granger_L.grangerspctrm = granger;
    
    M = mean(granger)
    r = mean(grangerflip.grangerspctrm)
    
    x = granger_L.freq(1:140);
    figure
    plot(x,M(1:140))
    hold on
    plot(x,r(1:140),':')
    xlabel('Frequency (Hz)')
    ylabel('Granger Causality')
    legend('Visual','Visual Flipped')
    title([subject{i}])
end

%% ASD

subject = {'0401', '0402','0404','0405','0406'};
granger_comb = [];   
for i=1:length(subject)
    cd(sprintf('D:\\ASD_Data\\%s\\visual\\granger',subject{i}'));
    load('granger_L.mat'); 
    load('granger_R.mat');
    granger = granger_R.grangerspctrm + granger_L.grangerspctrm;
    granger_R.grangerspctrm = granger;
    
    M = mean(granger)./2
    r = mean(grangerflip.grangerspctrm)
    
    x = granger_R.freq(35:140);
    figure
    plot(x,M(35:140))
    hold on
    plot(x,r(35:140),':')
    xlabel('Frequency (Hz)')
    ylabel('Granger Causality')
    legend('Visual','Visual Flipped')
    title([subject{i}])
    
    x = granger_R.freq(1:35);
    figure
    plot(x,M(1:35))
    hold on
    plot(x,r(1:35),':')
    xlabel('Frequency (Hz)')
    ylabel('Granger Causality')
    legend('Visual','Visual Flipped')
    title([subject{i}])
end