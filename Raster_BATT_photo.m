function Raster_BATT_photo(varargin)
%% AUTHOR    : Oscar X. Guerrero-Gutierrez
%% $DATE     : 02-Mar-2019 $
%% DEVELOPED : (R2015a)
%% FILENAME  : Raster_BATT_photo.m

%% Set defaults and load optional arguments
% argument parser
p = inputParser;
p.KeepUnmatched = true;
% parameter definition
% required parameters
p.addRequired('data');
% optional parametes
p.addParameter('slcn_seq', 0);
p.addParameter('minlag', -3);
p.addParameter('maxlag', 6);
p.addParameter('resolution', 0.05);
p.addParameter('rwd_period', 4);
p.addParameter('color_licks_out', [211,211,211]/255);
p.addParameter('color_rwd', [180,67,212]/255);  % "Lila"
p.addParameter('color_line', [211,211,211]/255);  % Grey
p.addParameter('colorStruct', 0);
p.addParameter('block_height', 0.5);
p.addParameter('spacing_bout_type', 0.3);
p.addParameter('spacing_slcn', 0.3);
p.addParameter('fig_title', []);
p.addParameter('raster_title', []);
% smooth parameters
p.addParameter('L', 5);

% read arguments
p.parse(varargin{:})
args = p.Results;
% TODO: check data structure
% dependent default values
if ~iscell(args.slcn_seq)
    args.slcn_seq = fieldnames(args.data.Solutions)';
end
slcn_seq = args.slcn_seq;

if ~isstruct(args.colorStruct)
    args.colorStruct = struct;
%     default_colors = {'b', 'g', 'y', 'r', 'k'};
    for i_slcn = 1:numel(slcn_seq)
        slcn = slcn_seq{i_slcn};
        args.colorStruct.(slcn) = 'k';
    end
end

[~, filename] = fileparts(args.data.Parameters.File);
if ~ischar(args.fig_title)
    args.fig_title = strcat(filename, ' - Raster');
end
if ~ischar(args.raster_title)
    args.raster_title = filename;
end
args.raster_title = strrep(args.raster_title, '_', '\_');

% plot parameters
plotArgs.LineWidth = 2;
plotArgs.LineStyle = '-';
% plotArgs.Color = P.color_licks_out;
tmp = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)];
otherArgs = reshape(tmp',[],1)';
% TODO: complete update
for arg_value = otherArgs
    [arg, value] = arg_value{:};
    plotArgs.(arg) = value;
end

args.lags = lagsfun(args.minlag, args.maxlag, args.resolution);

% aliases
data = args.data;
Solutions = data.Solutions;
trialsTable = data.Trials;
Rwd = data.Values.Rwd;
fig_title = args.fig_title;
raster_title = args.raster_title;
colors = args.colorStruct;
colors.licks_out = args.color_licks_out;
colors.rwd = args.color_rwd;
colors.line = args.color_line;
lags = args.lags;

has_photo = isfield(data.Values, 'photo');
if has_photo
    photo_time = data.Values.photo.('time');
    n_x_photo = max(cellfun(@length, trialsTable.dF));
    x_photo = linspace(args.minlag, args.maxlag, n_x_photo);
end

% function to generate raster values for each trial
x_values = @(x) repmat(x', 2, 1);
y_values = @(n, j, k) ([-(j+k) 0; -(j+k) args.block_height] * ones(2, n));

%% main
figure('NumberTitle', 'off', 'Name', fig_title);

hist_photo = cell(numel(slcn_seq), 1);

% raster initialization
if has_photo
    subplot(5, 1, 1:3);
else
    subplot(4, 1, 1:3);
end
hold on
ylabel('Trials');
title(raster_title);

% iterate over solution names
k = 0;  % spacing counter
for i_slcn = 1:numel(slcn_seq)
    slcn = slcn_seq{i_slcn};
    % get table of all trials for `slcn`
    trials_slcn = Solutions.(slcn).trials;
    % sort trials by bout type (incomplete or complete)
    for bout_type = {'incomplete' 'complete'}
        bout_type = char(bout_type);
        % filter trials of bout type
        il_trial_type = strcmp(trials_slcn.type, bout_type);
        i_trial_type = find(il_trial_type);
        for j = 1:numel(i_trial_type)
            i_trial = i_trial_type(j);
            % get values and init from solution trials table
            trial_values = trials_slcn.values{i_trial};
            init = trials_slcn.event(i_trial);
            
            % plot trial_values
            rel_trial_values = trial_values - init;
            n_values = length(trial_values);
            xx_trial = x_values(rel_trial_values);
            yy_trial = y_values(n_values, j, k);
            plotArgs.Color = colors.(slcn);

            plot(xx_trial, yy_trial, plotArgs);
            
            % plot licks off-trial
            % get trial window
            trial_licks = trials_slcn.licks{i_trial};
            il_offtrial = ~ismember(trial_licks, trial_values);
            licks_offtrial = trial_licks(il_offtrial);
            
            rel_licks = licks_offtrial - init;
            n_licks = length(licks_offtrial);
            xx_offtrial = x_values(rel_licks);
            yy_offtrial = y_values(n_licks, j, k);
            plotArgs.Color = colors.licks_out;
            
            plot(xx_offtrial, yy_offtrial, plotArgs);
            
            % plot accessible reward
            init_window = init + args.minlag;
            fin_window = init + args.maxlag;
            il_rwd = (init_window <= Rwd) & (Rwd <= fin_window);
            if any(il_rwd)
                trial_rwd = Rwd(il_rwd) - init;
                n_rwd = length(trial_rwd);
                xx_rwd = x_values(trial_rwd);
                yy_rwd = y_values(n_rwd, j, k);
                plotArgs.Color = colors.rwd;
                
                plot(xx_rwd, yy_rwd, plotArgs);
            end
            
            if has_photo
                photo_values = trials_slcn.dF{i_trial}';
                % pad trials with nan to make them stackable
                n_values = length(photo_values);
                n_diff = n_x_photo - n_values;
                if n_diff
                    % in case only one value is missing append a nan
                    if n_diff == 1
                        photo_values = [photo_values nan];
                    % in case first trial has less values
                    elseif trials_slcn.index(i_trial) == 1
                        photo_values = [nan(1, n_diff) photo_values];
                    % otherwise assume session ended before last trial
                    else
                        photo_values = [photo_values nan(1, n_diff)];
                    end
                end
                % build histogram data
                hist_photo{i_slcn} = [hist_photo{i_slcn}; photo_values];
                % plot photo line
                y = photo_values;
                rel_photo = (y - min(y)) / (max(y) - min(y));
                y = args.block_height - rel_photo - (j+k);
                plot(x_photo, y, 'Color', 'k');

            end  % if photo

            
        end  % for trial
        % in case we have an emtpy incomplete/complete trial vector
        if isempty(j)
            j = 0;
        end
        k = k + j + args.spacing_bout_type;
    end  % for bout_type
    % draw line between solutions
    if i_slcn ~= numel(slcn_seq)
        xx_line = [args.minlag; args.maxlag];
        y_line = k + args.spacing_slcn;
        yy_line = [-y_line; -y_line];
        plotArgs.Color = colors.line;
        plot(xx_line, yy_line, plotArgs);
    end
    k = k + args.spacing_slcn;
end  % for solution


% plot lick PSTH
if has_photo
    subplot(5, 1, 4);
else
    subplot(4, 1, 4);
end
hold on
ylabel('Licks/s');
for i_slcn = 1:numel(slcn_seq)
    slcn = slcn_seq{i_slcn};
    licks_all = data.Values.Licks;
    slcn_events = Solutions.(slcn).trials.event;
    [n,~,~,~] = PSTH2(licks_all, slcn_events, args.minlag, args.maxlag, args.resolution);
    sn = smoothPSTH(n, 3, args.lags);
    plot(lags, sn, 'Color', colors.(slcn));
end


% plot photo histogram
if has_photo
    % initialize plot
    subplot(5, 1, 5);
    hold on
    ylabel('dF/F');
    xlabel('Time');
    %
    for i_slcn = 1:numel(slcn_seq)
        slcn = slcn_seq{i_slcn};
        x = x_photo;
        y = nanmean(hist_photo{i_slcn}, 1);
        y = smoothPSTH(y, args.L);
        plot(x, y, 'Color', colors.(slcn));

    end
else
    warning('Photo data not found, will not be plotted');
end


end