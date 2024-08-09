clear all; close all;clc;

directory_path='C:\Users\user\Desktop\CommunicationSurprise_DataCode';%please make sure that your working directory is set to ...\CommunicationSurprise_DataCode
cd(directory_path); 

load('Participants-generated and Model-simulated Messages\behavioural and simulated data\percentage_data1.mat'); %import percentage data for the 4 types of the dataset( if you want to run these code for dataset 2, change data1 to data 2
% D is a 116x19 table and it includes 4 types of datas that are identified
% with group number
%group 1: participants data
%group 2: SMM
%group 3: SM
%group 4: MM

bf_mtype = readmatrix('Participants-generated and Model-simulated Messages\bayes factors\Bf_mtype1_surprise_models.csv') %import BF for mtype dataset (change mtype1 to mtype2 if you want to run analysis for the second dataset
bf_ttype = readmatrix('Participants-generated and Model-simulated Messages\bayes factors\Bf_ttype1_surprise_models.csv') %import BF for trial type
bf_mfeature = readmatrix('Participants-generated and Model-simulated Messages\bayes factors\Bf_mfeature1_surprise_models.csv') %import BF for the message profile (features)

group_id=1; % by changing this group identifier you can choose which dataset you plot. For exmaple if you want to plot behavioural data you choose 1 (as it is given in the exmaple now). Refer to group ids above to get the correct data

plot_mtype(D,group_id,bf_mtype) %1) plot mtype
plot_ttype(D,group_id,bf_ttype)  %2) plot ttype
plot_mfeature(D,group_id,bf_mfeature) %3) plot mfeature



%%%%%%%%%%%%%%%%%%%%% sub-functions  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 1) plot mtype
function plot_mtype(data,id,bf)

ind_group=find(data.group==id);
data=data(ind_group,:);

%caclulate SE
SEM(1,1) = std(data.E(:,1))/sqrt(length(data.E(:,1)));
SEM(1,2) = std(data.W(:,1))/sqrt(length(data.W(:,1)));
SEM(1,3) = std(data.P(:,1))/sqrt(length(data.P(:,1)));
%calculate mean
meanE=mean(data.E(:,1));meanW=mean(data.W(:,1));meanP=mean(data.P(:,1));

%plot
figure(6)
ax = gca();
x = [meanE,meanW,meanP]; x=round(x);
legends={'Enter-Exit';'Wiggly';'Pass-By'}
colors = {[224, 122, 95]/256;   %hot pink
    [61, 64, 91]/256;   %spring green
    [129, 178, 154]/256};  %dark orchid
y = [1:3];
b=bar(y,x,'FaceColor','flat','BarWidth', 0.75);
b.CData(1,:) = colors{1,:};
b.CData(2,:) = colors{2,:};
b.CData(3,:) = colors{3,:};
ylim([0 100]);
% Calculate center of each bar
set(b,'LineWidth',0.01,'EdgeColor',[ 1     1     1 ]);
ylabel('Percentage %');

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:10:100, ...
    'xticklabel',legends)
hold on
errorbar(x, SEM, 'k', 'linestyle', 'none','CapSize',4,'Color', [.3 .3 .3]);
ctr = bsxfun(@plus, b(1).XData, [b(1).XOffset]');

for i=2:length(bf(:,id))
    if bf(i,id)<1 & x(i-1)~=0
        plot(ctr(i-1), x(i-1)+6,  '*','MarkerSize',6,'Color', [.3 .3 .3]);
    end
end

end

%%%%%%%%%%%%% 2) plot ttype
function plot_ttype(data,id,bf)
%get ttype data
ind_group=find(data.group==id);
data=data(ind_group,:);

e(:,1)=data.E_easy(:,1);
e(:,2)=data.W_easy(:,1);
e(:,3)=data.P_easy(:,1);
h(:,1)=data.E_hard(:,1);
h(:,2)=data.W_hard(:,1);
h(:,3)=data.P_hard(:,1);

%caclulate SE
SEM(1,:) = nanstd(e)/sqrt(length(e));
SEM(2,:) = nanstd(h)/sqrt(length(h));
%calculate mean
s(1,:)=nanmean(e);
s(2,:)=nanmean(h);

%plot
b=bar(s,'FaceColor','flat','BarWidth', 0.9);
colors = {[224, 122, 95]/256;   %hot pink
    [61, 64, 91]/256;   %spring green
    [129, 178, 154]/256};  %dark orchidset(b,{'FaceColor'},colors);
set(b,'LineWidth',0.01,'EdgeColor',[ 1     1     1]);
set(b,{'FaceColor'},colors);

legends={'Enter-Exit';'Wiggly';'Pass-By'}% model data

name = {'Indirect trial','Direct trial'};
ylim([0 100]);
ylabel('Percentage');
hold on

% Calculate center of each bar
[ngroups, nbars] = size(s);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, s(:,i), SEM(:,i), 'Color', [.3 .3 .3],'linestyle', 'none','CapSize',3);
end

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:10:100, ...
    'xticklabel',name)
hold on

ctr1 = bsxfun(@plus, b(1).XData, [b(1).XOffset]');
ctr2 = bsxfun(@plus, b(2).XData, [b(2).XOffset]');
ctr3 = bsxfun(@plus, b(3).XData, [b(3).XOffset]');

ctr=[ctr1,ctr2,ctr3];
ctr=sort(ctr);
s=s';s=s(:);

for i=2:length(bf(:,id))
    if bf(i,id)<1 & s(i-1)~=0
        plot(ctr(i-1), s(i-1)+6,  '*','MarkerSize',6,'Color', [.3 .3 .3]);
    end
end
legend(legends,'Orientation','vertical','TextColor',[.3 .3 .3]);
legend boxoff  ;

end
%%%%%%%%%%%%% 3) plot mfeature

function plot_mfeature(data,id,bf)
% get the mfeature data
ind_group=find(data.group==id);
data=data(ind_group,:);

E(:,1)=data.E_Forward(:,1);
E(:,2)=data.E_backward(:,1);
E(:,3)=data.E_left_right(:,1);

W(:,1)=data.W_forward(:,1);
W(:,2)=data.W_backward(:,1);
W(:,3)=data.W_left_right(:,1);

P(:,1)=data.P_forward(:,1);
P(:,2)=data.P_backward(:,1);
P(:,3)=data.P_left_right(:,1);

%calculate the SE
SEM(1,:) = nanstd(E)/sqrt(length(E));
SEM(2,:) = nanstd(W)/sqrt(length(W));
SEM(3,:) = nanstd(P)/sqrt(length(P));

%calculate the mean
s(1,:)=nanmean(E);
s(2,:)=nanmean(W);
s(3,:)=nanmean(P);

%plot
b=bar(s,'FaceColor','flat','BarWidth', 0.9);
colors = {[130, 112, 129]/256;   %hot pink
    [222, 186, 111]/256;   %spring green
    [173, 106, 108]/256};  %darkrk orchid
set(b,{'FaceColor'},colors);
set(b,'LineWidth',0.01,'EdgeColor',[ 1     1     1]);
legends={'Forward';'Backward';'Left/Right'}

name = {'Enter-Exit';'Wiggly';'Pass-By'};
ylabel('Percentage');
hold on

% Calculate center of each bar
[ngroups, nbars] = size(s);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, s(:,i), SEM(:,i), 'Color', [.3 .3 .3],'linestyle', 'none','CapSize',3);
end
ylim([0 40]);
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'YTick'       , 0:10:100, ...
    'xticklabel',name)
hold on


ctr1 = bsxfun(@plus, b(1).XData, [b(1).XOffset]');
ctr2 = bsxfun(@plus, b(2).XData, [b(2).XOffset]');
ctr3 = bsxfun(@plus, b(3).XData, [b(3).XOffset]');
ctr=[ctr1,ctr2,ctr3];
ctr=sort(ctr);
s=s';s=s(:);

for i=2:length(bf(:,id))
    if bf(i,id)<1 & s(i-1)~=0
        plot(ctr(i-1), s(i-1)+15,  '*','MarkerSize',6,'Color', [.3 .3 .3]);
    end
end
legend(legends,'Orientation','vertical','TextColor',[.3 .3 .3]); % indicate confidence sender b2
legend boxoff  ;


end