
clear all; close all; clc;

directory_path = 'C:\Users\user\Desktop\CommunicationSurprise_DataCode'; % Set the correct working directory
cd(directory_path); 

load('Participants-generated and Model-simulated Messages\behavioural and simulated data\percentage_data2.mat'); % Import data
bf_mtype = readmatrix('Participants-generated and Model-simulated Messages\bayes factors\Bf_mtype1_surprise_models.csv'); % Import BF for mtype
bf_ttype = readmatrix('Participants-generated and Model-simulated Messages\bayes factors\Bf_ttype1_surprise_models.csv'); % Import BF for trial type
bf_mfeature = readmatrix('Participants-generated and Model-simulated Messages\bayes factors\Bf_mfeature1_surprise_models.csv'); % Import BF for features

group_id = 1; % Set group identifier

% Uncomment as needed to plot different features
plot_mtype(D, group_id, bf_mtype); % Plot mtype
plot_ttype(D, group_id, bf_ttype); % Plot ttype
plot_mfeature(D, group_id, bf_mfeature); % Plot mfeature

%%%%%%%%%%%%%%%%%%%%% Sub-Functions %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% 1) plot mtype
function plot_mtype(data, id, bf)
ind_group = find(data.group == id);
data = data(ind_group, :);

% Calculate SE and means
SEM = [std(data.E(:,1)), std(data.W(:,1)), std(data.P(:,1))] ./ sqrt(size(data.E, 1));
means = [mean(data.E(:,1)), mean(data.W(:,1)), mean(data.P(:,1))];
%means = [mean(nonzeros(data.E(:,1))), mean(nonzeros(data.W(:,1))), mean(nonzeros(data.P(:,1)))];

% Plot bar graph
figure;
x = means;
legends = {'EE', 'W', 'PB'};
colors = {[224, 122, 95]/256, [61, 64, 91]/256, [129, 178, 154]/256}; % Base colors
light_colors = cellfun(@(c) c + 0.5 * (1 - c), colors, 'UniformOutput', false); % Lighter colors
y = 1:3;
b = bar(y, x, 'FaceColor', 'flat', 'BarWidth', 0.75);
for i = 1:length(colors)
    b.CData(i,:) = colors{i};
end
ylim([0 100]);
set(b, 'LineWidth', 0.01, 'EdgeColor', [1 1 1], 'FaceAlpha', 0.7); % Adjust transparency
ylabel('%');
set(gca, 'Box', 'off', 'TickDir', 'out', 'XMinorTick', 'off', ...
    'YMinorTick', 'on', 'YGrid', 'off', 'XColor', [.3 .3 .3], ...
    'YColor', [.3 .3 .3], 'YTick', 0:20:100, 'xticklabel', legends);
hold on;


% Overlay individual data points with jitter
x_offset = 0.4; % Adjust horizontal spread of data points
for i = 1:length(colors)
    % Select corresponding data field dynamically
    if i == 1
        field_data = data.E(:, 1); % Enter-Exit
    elseif i == 2
        field_data = data.W(:, 1); % Wiggly
    elseif i == 3
        field_data = data.P(:, 1); % Pass-By
    end
    scatter(i + x_offset * (rand(size(field_data, 1), 1) - 0.5), field_data, 10, 'filled', ...
        'MarkerFaceColor', light_colors{i}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
end

% Significance stars
ctr = y;
for i = 2:length(bf(:,id))
    if bf(i,id) < 1 && x(i-1) ~= 0
        plot(ctr(i-1), x(i-1)+6, '*', 'MarkerSize', 6, 'Color', [.3 .3 .3]);
    end
end

hold on;

% Error bars
errorbar(y, x, SEM, 'k', 'linestyle', 'none', 'CapSize', 4, 'Color', [.3 .3 .3]);
end

%%%%%%%%%%%%% 2) plot ttype
function plot_ttype(data, id, bf)
ind_group = find(data.group == id);
data = data(ind_group, :);

% Prepare data for plotting
e = [data.E_easy(:,1), data.W_easy(:,1), data.P_easy(:,1)];
h = [data.E_hard(:,1), data.W_hard(:,1), data.P_hard(:,1)];
SEM = [nanstd(e)/sqrt(size(e, 1)); nanstd(h)/sqrt(size(h, 1))];
means = [nanmean(e); nanmean(h)];

% Plot bar graph
figure;
b = bar(means, 'FaceColor', 'flat', 'BarWidth', 0.9,'FaceAlpha', 0.7,'EdgeColor', [1 1 1]);
colors = {[224, 122, 95]/256, [61, 64, 91]/256, [129, 178, 154]/256}; % Base colors
light_colors = cellfun(@(c) c + 0.5 * (1 - c), colors, 'UniformOutput', false); % Lighter colors
for i = 1:length(colors)
    b(i).FaceColor = colors{i};
end
ylim([0 100]);
ylabel('%');

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:20:100, ...
    'xticklabel', {'Indirect trial', 'Direct trial'});
hold on;

% Error bars
[ngroups, nbars] = size(means);
groupwidth = min(0.8, nbars/(nbars + 1.5));


% Overlay individual data points with jitter
x_offset = 0.17; % Adjust horizontal spread of data points
for i = 1:ngroups
    for j = 1:nbars
        % Dynamically select data fields for individual points
        if i == 1
            field_data = e(:, j); % Easy trials
        elseif i == 2
            field_data = h(:, j); % Hard trials
        end
        x_jitter = (i - groupwidth/2 + (2*j-1) * groupwidth / (2*nbars)) + x_offset * (rand(size(field_data, 1), 1) - 0.5);
        scatter(x_jitter, field_data, 10, 'filled', ...
            'MarkerFaceColor', light_colors{j}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3); % Match colors with bars
    end
end
hold on
for j = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*j-1) * groupwidth / (2*nbars);
    errorbar(x, means(:,j), SEM(:,j), 'Color', [.3 .3 .3], 'linestyle', 'none', 'CapSize', 3);
end
% Significance stars
ctr = sort(b(1).XData);
for i = 2:length(bf(:,id))
    if bf(i,id) < 1
        plot(ctr(i-1), means(i-1)+6, '*', 'MarkerSize', 6, 'Color', [.3 .3 .3]);
    end
end
legend({'Enter-Exit', 'Wiggly', 'Pass-By'}, 'Location', 'northeast', 'TextColor', [.3 .3 .3]);
legend boxoff;
end

%%%%%%%%%%%%% 3) plot mfeature
function plot_mfeature(data, id, bf)
ind_group = find(data.group == id);
data = data(ind_group, :);

% Prepare data for plotting
E = [data.E_Forward(:,1), data.E_backward(:,1), data.E_left_right(:,1)];
W = [data.W_forward(:,1), data.W_backward(:,1), data.W_left_right(:,1)];
P = [data.P_forward(:,1), data.P_backward(:,1), data.P_left_right(:,1)];
SEM = [nanstd(E)/sqrt(size(E, 1)); nanstd(W)/sqrt(size(W, 1)); nanstd(P)/sqrt(size(P, 1))];
means = [nanmean(E); nanmean(W); nanmean(P)];

% Plot bar graph
figure;
b = bar(means, 'FaceColor', 'flat', 'BarWidth', 0.9,'FaceAlpha', 0.7,'EdgeColor', [1 1 1]);
colors = {[130, 112, 129]/256, [222, 186, 111]/256, [173, 106, 108]/256};
light_colors = cellfun(@(c) c + 0.5 * (1 - c), colors, 'UniformOutput', false);
for i = 1:length(colors)
    b(i).FaceColor = colors{i};
end
ylim([0 60]);
ylabel('%');
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:10:60, ...
    'xticklabel', {'EE', 'W', 'PB'});
hold on;

% Error bars
[ngroups, nbars] = size(means);
groupwidth = min(0.8, nbars/(nbars + 1.5));


% Overlay individual data points with jitter
x_offset = 0.14;
for i = 1:ngroups
    for j = 1:nbars
        % Dynamically select data fields
        if i == 1
            field_data = E(:, j);
        elseif i == 2
            field_data = W(:, j);
        elseif i == 3
            field_data = P(:, j);
        end
        x_jitter = (i - groupwidth/2 + (2*j-1) * groupwidth / (2*nbars)) + x_offset * (rand(size(field_data, 1), 1) - 0.5);
        scatter(x_jitter, field_data, 10, 'filled', ...
            'MarkerFaceColor', light_colors{j}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3);
    end
end
hold on
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, means(:,i), SEM(:,i), 'Color', [.3 .3 .3], 'linestyle', 'none', 'CapSize', 3);
end

end
