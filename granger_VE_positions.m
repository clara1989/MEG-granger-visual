%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to show the positon of the atlas points used for GC analysis
%
% Written by Robert Seymour February 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get the Atlas Points of Interest on the HCP atlas

% Get ROIs based on L/R hemisphere

load('D:\scripts\MEG-granger-visual\atlas_MSMAll_4k.mat');

subject = {'RS'}; i = 1;

%% Load variables required for source analysis
load(sprintf('D:\\pilot\\%s\\visual\\data_clean_noICA.mat',subject{i}))
load(sprintf('D:\\pilot\\%s\\visual\\sourceloc\\mri_realigned.mat',subject{i}))
load(sprintf('D:\\pilot\\%s\\visual\\sourceloc\\sens.mat',subject{i}))
load(sprintf('D:\\pilot\\%s\\visual\\sourceloc\\seg.mat',subject{i}))
data_filtered = data_clean_noICA;

%% Set the current directory
cd(sprintf('D:\\pilot\\%s\\visual\\granger',subject{i}))

sourcespace = ft_read_headshape({['D:\megconnectome-3.0\template\Conte69.L.inflated.4k_fs_LR.surf.gii']});%,['D:\megconnectome-3.0\template\Conte69.R.inflated.4k_fs_LR.surf.gii']});
%% Set bad channel list - can change to specific channels if necessary
chans_included = {'MEG', '-MEG0322', '-MEG2542','-MEG0111','-MEG0532'};

%% Load 3D 4k Cortical Mesh for L/R hemisphere & Concatenate

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

%% Plot the positions of the VE for each visual ROI
pntV1 = find(atlas.parcellation==2);%pntV1 = pntV1-4002; 
pntV1 = pntV1(5:49); 
pntV4 = find(atlas.parcellation==7);%pntV4 = pntV4-4002; 
pntV4 = pntV4(1:13); 

figure; hold on;
dataV1 = ft_plot_mesh(sourcemodel_virt.pos(pntV1,:),'vertexcolor','y');
dataV4 = ft_plot_mesh(sourcemodel_virt.pos(pntV4,:),'vertexcolor','g');
legend([dataV1 dataV4],'V1','V4'); hold on;
ft_plot_mesh(sourcespace,'facecolor','k','edgecolor','none','facealpha',0.1); camlight;
set(gcf,'color','w'); view(-30, 0); drawnow;

