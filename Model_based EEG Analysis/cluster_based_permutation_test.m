%%%% analyzing the  regression weights by testing these weights against 0 across subjects for all electrodes and time-points.
clear all; close all; clc;

%% 1) set working directory,load the data, and initiate fieldtrip
directory_path='C:\Users\user\Desktop\CommunicationSurprise_DataCode';%please make sure that your working directory is set to ...\CommunicationSurprise_DataCode
cd(directory_path); 

addpath('D:\TCG_project\EEG-analysis\fieldtrip-20220104'); %% add fieldtrip to the path in order to run this code
addpath(genpath('D:\TCG_project\EEG-analysis\eeglab\eeglab2021.1\functions')); %% add eeglab to the path. We only need eeglab to use the readlocs function)
ft_defaults; %initiate Fieldtrip

loc_data_path='Model_based EEG Analysis\data\';
filename.loc = '128 EQ A on MNI Colin27 with ABCD labels for BESA';
filename.elec = strcat(loc_data_path,filename.loc,'.sfp');
elec = readlocs(filename.elec); %read the electrode locations from the template
load('Model_based EEG Analysis\data\regression weights.mat'); %load the regression weights
load('Model_based EEG Analysis\data\meta_data\time.mat'); %load the vector of time of the observational section
load('Model_based EEG Analysis\data\meta_data\labels.mat');  %laod the labels

%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%
%% 2) create the data structures that can work with fieldtrip for regression weigths data and mock data
for subindx=[1:5,7:10,12:23,25:32]

    beta=b(subindx,:);
    beta=vertcat(beta{:});

    % beta time series
    avg{subindx}.avg=beta;
    avg{subindx}.label=label;
    avg{subindx}.time=time;

    % zero time series (mock data contains zeroes but has the same
    % structure as the regression weights data
    avg_mock{subindx}=avg{subindx};
    avg_mock{subindx}.avg=zeros(128,307);
end

%%  3) get the layout from the location file
cfg = [];
cfg.layout = filename.elec;
layout = ft_prepare_layout(cfg, avg{1});
layout.pos(:,1)=layout.pos(:,1)*1.13;
layout.pos(:,2)=layout.pos(:,2)*1.3;
layout.pos(:,2)=layout.pos(:,2)+0.07;

%% 4) delete the empty cells
avg=avg(~cellfun('isempty',avg));
avg_mock=avg_mock(~cellfun('isempty',avg_mock));

%% 5) find the grand avarage for both regression and mock data
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
grandavg      = ft_timelockgrandaverage(cfg, avg{:});
grandavg_mock      = ft_timelockgrandaverage(cfg, avg_mock{:});

%% 6) get statistics for grandavarage (ft_timelockstatistics)
% now create neighbors
cfg_neighb        = [];
cfg_neighb.method = 'distance';
cfg_neighb.layout = layout;
neighbours        = ft_prepare_neighbours(cfg_neighb);

% create new cfg
cfg         = [];
cfg.channel = {'EEG'};
cfg.latency = 'all';
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';%'indepsamplesT'
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.01;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.neighbours       = neighbours;  % same as defined for the between-trials experiment
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
Nsubj  = 29;
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];
cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;
[stat] = ft_timelockstatistics(cfg, avg{:}, avg_mock{:});
stat.posclusters(1);
stat.posclusterslabelmat

%% 6) run ft_math
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_betavszero    = ft_math(cfg,grandavg, grandavg_mock);


%% 7) define parameters for plotting
timestep      =   0.2; %(in seconds)
sampling_rate = 256;
sample_count  = length(stat.time);
j = [0:timestep:1.2];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:51:sample_count];  % temporal endpoints in EEG samples
% get relevant positvie values
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(stat.posclusterslabelmat, pos_clust);
% get negative values
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_clust = find(neg_cluster_pvals < 0.025);
neg       = ismember(stat.negclusterslabelmat, neg_clust);
[i1,i2] = match_str(GA_betavszero.label, stat.label);

%% 8) plot
for k = 1:6
    cfg=[];
    cfg.figure     = subplot(4,4,k);
    cfg.xlim       = [round(j(k),3) round(j(k+1),3)];
    pos_int        = zeros(numel(GA_betavszero.label),1);
    pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
    cfg.highlight  = 'on';
    cfg.highlightsymbol  = {'.'};
    cfg.highlightcolor='w';
    cfg.highlightsize    = 8;
    cfg.highlightchannel = find(pos_int);
    cfg.comment    = 'xlim';
    cfg.commentpos = 'title';
    cfg.layout     = layout;
    cfg.figure     = 'gca';
    cfg.colormap='*RdBu';
    cfg.colorbar='yes';
    %cfg.zlim=[-0.6,0.6]
    ft_topoplotER(cfg, GA_betavszero);
    sign_ell{k}=find(pos_int);
end

%% 9) plot the significent electrode clusters
% get the significant channels and their corresponding regression data
front_chans_group=sign_ell(3:4);
front_chans_all=[front_chans_group{1};front_chans_group{2}]
front_chans_all=unique(front_chans_all);
central_chans_group=sign_ell(5:6);
central_chans_all=[central_chans_group{1};central_chans_group{2}]
central_chans_all=unique(central_chans_all);
common_frontal_ind=ismember(front_chans_all,central_chans_all);
front_chans_all(common_frontal_ind)=[];
roi_front_betas=GA_betavszero.avg(front_chans_all,:);
avg_front_betas=mean(roi_front_betas);

% labels
labels_frontal=labels(front_chans_all);
labels_central=labels(central_chans_all);
% smooth the data with factor 15
windowsize=15;
avg_smooth_frontal=smooth(avg_front_betas,windowsize);
%plot the frontal cluster regressors
p=plot(time,avg_smooth_frontal, 'Color', [0.3765    0.6667    0.7412]);
selected_area = m(3):m(5);
hold on
area(time(selected_area),avg_smooth_frontal(m(3):m(5))-0.01,'EdgeColor','none','FaceColor', [ 0.9686    0.8157    0.8157],'LineStyle',':');
yli=[-0.4,0.9];yti=[-0.4:0.2:0.9];
xli=[0,1.2];xti=j;
set(gca, 'YTick',yti) ;
ylim(yli);
set(gca, 'XTick',xti) ;
xlim(xli);
box off

% plot the central-frontal cluster electrodes
roi_central_betas=GA_betavszero.avg(central_chans_all,:);
avg_central_betas=mean(roi_central_betas);
avg_smooth_central=smooth(avg_central_betas,windowsize);
p=plot(time,avg_smooth_central, 'Color', [0.3765    0.6667    0.7412]);
selected_area = m(5):m(7);
hold on
area(time(selected_area),avg_smooth_central(m(5):m(7))-0.01,'EdgeColor','none','FaceColor', [ 0.9686    0.8157    0.8157],'LineStyle',':');
yli=[-0.4,0.9];yti=[-0.4:0.2:0.9];
xli=[0,1.2];xti=j;
set(gca, 'YTick',yti) ;
ylim(yli);
set(gca, 'XTick',xti) ;
xlim(xli);
box off

