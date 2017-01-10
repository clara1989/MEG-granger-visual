%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script uses the GC spectra calculated using the
% granger_visual_ASD.m script, sorts into feedforward vs feedback
% connections (defined using macaque data), concatenates and plots 
% ASD vs control
%
% Statistical analysis to come...
%
% Written by Robert Seymour January 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load in data

subject = {'0401','0402','0403','0404','0405','0406','0407','0409','0413','0414'};  
feedforward_ASD = [];
feedback_ASD = [];

for i=1:length(subject)
    cd(sprintf('D:\\ASD_Data\\%s\\visual\\granger',subject{i}'));
    load('granger_L.mat'); 
    load('granger_R.mat');
    granger = (granger_L.grangerspctrm + granger_R.grangerspctrm)./2;
    
    ff = [1,3,5,7,9,11,13,15,17,19,21,23,26,27,30];
    fb = [2,4,6,8,10,12,14,16,18,20,22,24,25,28,29];
    
    for j = 1:length(ff)
        feedforward_ASD = vertcat(feedforward_ASD,granger(ff(j),:));
        
    end
for k = 1:length(ff)
        feedback_ASD = vertcat(feedback_ASD,granger(fb(k),:));
    end
end


%% Control Data

subject = {'1401','1402','1403'};
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
        feedforward = vertcat(feedforward,granger(ff(j),:));
        
    end
for k = 1:length(ff)
        feedback = vertcat(feedback,granger(fb(k),:));
    end
end

%% CD to group folder
cd('D:\ASD_Data\Group\GC')

%% Plot feedforward_ASD vs feedback_ASD

mean_feedforward = mean(feedforward);
mean_feedback = mean(feedback);
mean_feedforward_ASD = mean(feedforward_ASD);
mean_feedback_ASD = mean(feedback_ASD);

% Separate Plots
x = granger_R.freq(1:140);
figure
subplot(2,1,1)
plot(x,mean_feedforward)
hold on
plot(x,mean_feedback_ASD)
xlabel('Frequency (Hz)')
ylabel('Granger Causality')
legend('feedforward ASD','feedback ASD')
title('ASD')
ylim([0 0.035]);
subplot(2,1,2)
hold on 
plot(x,mean_feedback);
hold on 
plot(x,mean_feedforward_ASD)
ylim([0 0.025]);
xlabel('Frequency (Hz)')
ylabel('Granger Causality')
title('Control')
legend('feedforward control','feedback Control')
saveas(gcf,'ff_fb_separate_plot.png')   

% Plots all together
x = granger_R.freq(1:140);
figure
plot(x,mean_feedforward_ASD)
hold on
plot(x,mean_feedback_ASD)
hold on 
plot(x,mean_feedforward);
hold on 
plot(x,mean_feedback);
xlabel('Frequency (Hz)')
ylabel('Granger Causality')
legend('feedforward ASD','feedback ASD','feedforward control','feedback Control')%,'GC flipped')
saveas(gcf,'ff_fb_plot_together.png')   

% Write files to CSV for external plotting
csvwrite('ASD_fb.csv',feedback_ASD); csvwrite('ASD_ff.csv',feedforward_ASD); 
csvwrite('control_fb.csv',feedback); csvwrite('control_ff.csv',feedforward); 

%% Error Bars

% CI = zeros(2,140);
% SEM = std(feedback)/sqrt(length(feedback));               % Standard Error
% ts = tinv([0.025  0.975],length(feedback)-1);      % T-Score
% CI(1,:) = mean(feedback) + rot90(ts(1).*SEM(:));        
% CI(2,:) = mean(feedback) + rot90(ts(2).*SEM(:));    % Confidence Intervals
% CI = CI./10;
% 
% figure
% shadedErrorBar(x,mean_feedback,CI(:,1:140))
% legend('95% Confidence Interval','Visual')
% xlabel('Frequency (Hz)')
% ylabel('Granger Causality')

%% DAI

DAI1 = zeros(size(feedback));

for df = 1:140
    DAI1(:,df) = ((feedforward(:,df) - feedback(:,df)) ./ (feedforward(:,df) + feedback(:,df)));
end

DAI2 = zeros(size(feedback_ASD));

for df = 1:140
    DAI2(:,df) = ((feedforward_ASD(:,df) - feedback_ASD(:,df)) ./ (feedforward_ASD(:,df) + feedback_ASD(:,df)));
end

x = [1:1:140];
figure
plot(x,mean(DAI1))
hold on
plot(x,mean(DAI2))



