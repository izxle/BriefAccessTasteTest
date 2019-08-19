function BATT_Descriptors(varargin)
%% Parameters
pArgs = inputParser;
pArgs.addRequired('dataTable');
pArgs.addRequired('slcn_seq');
%
pArgs.addParameter('color_slcn', 0);
pArgs.addParameter('color_group', 0);
pArgs.addParameter('marker_size', 110);
pArgs.addParameter('marker_seq', {'o'});
pArgs.addParameter('photo_x_time', [])
pArgs.addParameter('line_style_seq', {'-', '--', ':', '-.'});
pArgs.addParameter('photo_Zscore_y_lim', [-1.5 2.5])
pArgs.addParameter('photo_licks_y_lim', [0 14])
%
pArgs.addParameter('L', 5)
%
pArgs.addParameter('minlag', -3);
pArgs.addParameter('maxlag', 6);
pArgs.addParameter('resolution', 0.05);
%
pArgs.addParameter('bar_x_label', []);
%
pArgs.addParameter('reference_group', []);
%
pArgs.parse(varargin{:});
args = pArgs.Results;
% dependant default arguments
slcn_seq = args.slcn_seq;
if ~isstruct(args.color_slcn)
    args.color_slcn = struct;
    for i = 1:numel(args.slcn_seq)
        slcn = slcn_seq{i};
        args.color_slcn.(slcn) = 'k';
    end
end

if isempty(args.reference_group)
    % get last group as the reference group
    args.groups = unique(args.dataTable.Group);
    args.reference_group = args.groups{end};
end

if isempty(args.bar_x_label)
    args.bar_xlabel = slcn_seq;
end

% TODO: check marker_seq

lags = lagsfun(args.minlag, args.maxlag, args.resolution);

SEM = @(x)(nanstd(x)/sqrt(length(x)- sum(isnan(x))));

%% Alias
dataTable = args.dataTable;
color_group = args.color_group;
color_slcn = args.color_slcn;
marker_size = args.marker_size;
marker_seq = args.marker_seq;
line_style_seq = args.line_style_seq;
photo_x_time = args.photo_x_time;
photo_licks_y_lim = args.photo_licks_y_lim;
photo_Zscore_y_lim = args.photo_Zscore_y_lim;
reference_group = args.reference_group;

bar_x_label = args.bar_x_label;

%% Initialize variables
resTable = table;
%
PSTH = struct;
AUC = struct;
Correlations = struct;
%
groups = unique(dataTable.Group);
n_groups = length(groups);
n_slcns = length(slcn_seq);
windows = {'W1', 'W2', 'W3'}';
n_windows = length(fieldnames(dataTable.Data{1,1}.Correlations));
corr_info = {};

%% Re-group trials per subject
% regroup trials in resTable and add Subject and Solution labels
for i_group = 1:n_groups
    group = groups{i_group};
    il_table_group = strcmp(dataTable.Group, group);
    i_table_group = find(il_table_group);
    
    subjects_group = unique(dataTable.Subject(il_table_group));
    
    % initialize PSTH and AUC struct for group
    PSTH.(group) = struct;
    AUC.(group) = struct;
    
    for i = 1:n_slcns
        slcn = slcn_seq{i};
        PSTH.(group).(slcn).values = [];
        
        for i_window = 1:n_windows
            window = windows{i_window};

            AUC.(group).(window).(slcn).Zscore.values = [];
            AUC.(group).(window).(slcn).licks.values = [];
            AUC.(group).(window).(slcn).Z_max.values = [];
            AUC.(group).(window).(slcn).Z_peak_lat.values = [];
            AUC.(group).(window).line.x_licks = [];
            AUC.(group).(window).line.y_Zscore = [];
            AUC.(group).(window).line.y_Peak = [];
            AUC.(group).(window).line.y_Lat = [];
            Correlations.(group).(window).values = [];
            Correlations.(group).(window).fitline.Zs = [];
            Correlations.(group).(window).fitline.peak = [];
            Correlations.(group).(window).fitline.lat = [];
            Correlations.(group).(window).values_mean = [];
            Correlations.(group).(window).fitline.Zs_mean = [];
            Correlations.(group).(window).fitline.peak_mean = [];
            Correlations.(group).(window).fitline.lat_mean = [];
            Correlations.(group).(window).Results = [];
            Correlations.(group).(window).Results_mean = [];
            Correlations.Text.(window) = {};
       
            for i_g = 1:numel(subjects_group)
                subject = subjects_group{i_g};

                for i_window = 1:n_windows
                    window = windows{i_window};

                    AUC.(group).(window).subjects.(subject).(slcn).values = [];
                    AUC.(group).(window).subjects.(subject).(slcn).Zscore.values = [];
                    AUC.(group).(window).subjects.(subject).(slcn).licks.values = [];
                    AUC.(group).(window).subjects.(subject).(slcn).Z_max.values = [];
                    AUC.(group).(window).subjects.(subject).(slcn).Z_peak_lat.values = [];

                end
            end
        end
    end
    
    % iterate over files in each group
    for i = i_table_group'
        data = dataTable.Data{i};
        Solutions = data.Solutions;
        subject = dataTable.Subject{i};
        % iterate over solutions in each file
        for i_slcn = 1:n_slcns
            slcn = slcn_seq{i_slcn};
            
            trialsTable = Solutions.(slcn).trials;
            n_trials = size(trialsTable, 1);
            
            % add columns for group, subject and solution
            col_group = repmat({group}, n_trials, 1);
            trialsTable.Group = col_group;
            col_subject = repmat({subject}, n_trials, 1);
            trialsTable.Subject = col_subject;        
            col_solution = repmat({slcn}, n_trials, 1);
            trialsTable.Solution = col_solution;
            % add PSTH
            cell_licks = trialsTable.licks;
            licks = cat(1, cell_licks{:});
            events = trialsTable.event;
            
            [n, ~, ~, ~] = PSTH2(licks, events, args.minlag, args.maxlag, args.resolution);
            sn = smoothPSTH(n, 3, lags);
            PSTH.(group).(slcn).values = [PSTH.(group).(slcn).values; sn];
            
            for i_window = 1:n_windows
                window = windows{i_window};
                
                slcn_Z_AUC = Solutions.(slcn).PSTH.Zscore.(window).AUC_mean;
                slcn_licks_AUC = Solutions.(slcn).PSTH.licks.(window).AUC;
                slcn_Z_peak = Solutions.(slcn).PSTH.Zscore.(window).peak_mean;
                slcn_Z_latency = Solutions.(slcn).PSTH.Zscore.(window).latency_mean;
                
                AUC.(group).(window).(slcn).Zscore.values(end + 1) = slcn_Z_AUC;
                AUC.(group).(window).(slcn).licks.values(end + 1) = slcn_licks_AUC;
                AUC.(group).(window).(slcn).Z_max.values(end + 1) = slcn_Z_peak;
                AUC.(group).(window).(slcn).Z_peak_lat.values(end + 1) = slcn_Z_latency;
                  % save per subjects
                AUC.(group).(window).subjects.(subject).(slcn).Zscore.values(end + 1) = slcn_Z_AUC;
                AUC.(group).(window).subjects.(subject).(slcn).licks.values(end + 1) = slcn_licks_AUC;
                AUC.(group).(window).subjects.(subject).(slcn).Z_max.values(end + 1) = slcn_Z_peak;
                AUC.(group).(window).subjects.(subject).(slcn).Z_peak_lat.values(end + 1) = slcn_Z_latency;
                  % save regroup data for Pearson calculations(later)
                n_corr = size(Correlations.(group).(window).values, 1);
                Correlations.(group).(window).values(n_corr + 1, 1) = slcn_licks_AUC;
                Correlations.(group).(window).values(n_corr + 1, 2) = slcn_Z_AUC;
                Correlations.(group).(window).values(n_corr + 1, 3) = slcn_Z_peak;
                Correlations.(group).(window).values(n_corr + 1, 4) = slcn_Z_latency;
        
                % append trials to resTable
                resTable = [resTable; trialsTable];
                
            end
        end  % for solutions
    end  % for sessions in group
end  % for groups

% Pearson correlationships 
for i_group = 1:n_groups
    group = groups{i_group};
    
    for i_window = 1:n_windows
        window = windows{i_window};
        [rho, pval] = corr(Correlations.(group).(window).values);
        Correlations.(group).(window).Results.AUC_corr = rho(1,2);
        Correlations.(group).(window).Results.AUC_p = pval(1,2);
        Correlations.(group).(window).Text.r_AUC = strcat('r = ', num2str(rho(1,2)), ', p =', num2str(pval(1,2)));
        
        Correlations.(group).(window).Results.peak_corr = rho(1,3);
        Correlations.(group).(window).Results.peak_p = pval(1,3);
        Correlations.(group).(window).Text.r_peak = strcat('r = ', num2str(rho(1,3)), ', p =', num2str(pval(1,3)));
        
        Correlations.(group).(window).Results.latency_corr = rho(1,4);
        Correlations.(group).(window).Results.latency_p = pval(1,4);
        Correlations.(group).(window).Text.r_latency = strcat('r = ', num2str(rho(1,4)), ', p =', num2str(pval(1,4)));
        
        % Fit Line for Pearson correlation
        % get data points
        x = Correlations.(group).(window).values(:,1);
        y_Zs = Correlations.(group).(window).values(:,2);
        y_peak = Correlations.(group).(window).values(:,3);
        y_lat = Correlations.(group).(window).values(:,4);
        % Save beta coefficients for line
        Correlations.(group).(window).fitline.Zs = polyfit(x, y_Zs, 1);
        Correlations.(group).(window).fitline.peak = polyfit(x, y_peak, 1);
        Correlations.(group).(window).fitline.lat = polyfit(x, y_lat, 1);
        
    end % windows
 
end % groups 

% legends for correlations
for i_window = 1:n_windows
    window = windows{i_window};
    
    for i_group = 1:n_groups
        group = groups{i_group};
        
        Correlations.Text.(window).Text_AUC{i_group} =  Correlations.(group).(window).Text.r_AUC;
        Correlations.Text.(window).Text_peak{i_group} =  Correlations.(group).(window).Text.r_peak; 
        Correlations.Text.(window).Text_latency{i_group} =  Correlations.(group).(window).Text.r_latency;
        
    end
    
end

%% treat data

% initialize
data_lick_rate = struct;
data_bout_size = struct;
data_dF = struct;
data_Zscore = struct;

% pad trials with nan to make them stackable
% ignore any trial with missing data (incomplete baseline)
n_max = max(cellfun(@length, resTable.dF));
il_short = cellfun(@(x)(length(x) < n_max), resTable.dF);
for i = find(il_short)'
    n_i = length(resTable.dF{i});
    if n_max - n_i == 1
        vals = resTable.dF{i};
        resTable.dF{i} = [vals; nan];
        resTable.Zscore{i} = [vals; nan];
    elseif n_max - n_i == 2
        error('ERROR: debug here.')
    else
        resTable.dF{i} = nan(n_max, 1);
        resTable.Zscore{i} = nan(n_max, 1);
    end
end

% iterate over groups
for i_group = 1:n_groups
    group = groups{i_group};
    il_table_group = strcmp(resTable.Group, group);
    subjects_group = unique(resTable.Subject(il_table_group));
    plotArgs.LineStyle = line_style_seq{i_group};

    % iterate over solutions
    for i_slcn = 1:n_slcns
        slcn = slcn_seq{i_slcn};
        il_table_slcn = strcmp(resTable.Solution, slcn);
        il_table_group_slcn = il_table_group & il_table_slcn;
        % plot args
        barArgs.EdgeColor = color_slcn.(slcn);
        barArgs.FaceColor = 'none';
        scatterArgs.markerFaceColor = color_slcn.(slcn);
        scatterArgs.markerEdgeColor = color_slcn.(slcn);
        plotArgs.Color = color_slcn.(slcn);
        % get mean and SEM for
        % lick_rate
        lick_rate_values = resTable.lick_rate(il_table_group_slcn);
        data_lick_rate.(group).(slcn).mean = nanmean(lick_rate_values);
        data_lick_rate.(group).(slcn).SEM = SEM(lick_rate_values);
        data_lick_rate.(group).(slcn).barArgs = barArgs;
        % PSTH
        vals = PSTH.(group).(slcn).values;
        vals_mean = nanmean(vals, 1);
        vals_sem = nanstd(vals, 0, 1) / sqrt(size(vals, 1));
        data_lick_rate.(group).(slcn).PSTH.values = vals;
        data_lick_rate.(group).(slcn).PSTH.mean = vals_mean;
        data_lick_rate.(group).(slcn).PSTH.SEM = vals_sem;
        data_lick_rate.(group).(slcn).PSTH.plotArgs = plotArgs;
        
        % bout_size
        bout_size_values = resTable.bout_size(il_table_group_slcn);
        data_bout_size.(group).(slcn).mean = nanmean(bout_size_values);
        data_bout_size.(group).(slcn).SEM = SEM(bout_size_values);
        data_bout_size.(group).(slcn).barArgs = barArgs;
        
        % iterate over subjects
        for i_subject = 1:numel(subjects_group)
            subject = subjects_group{i_subject};
            
            il_table_subject = strcmp(resTable.Subject, subject);
            il_table_values = il_table_subject & il_table_slcn;
            
            % plot args
            scatterArgs.Marker = marker_seq{i_subject};
            
            % get mean values for each subject
            % lick rate
            lick_rate_vals = resTable.lick_rate(il_table_values);
            mean_val = nanmean(lick_rate_vals);
            data_lick_rate.(group).(slcn).Subjects.(subject).mean = mean_val;
            data_lick_rate.(group).(slcn).Subjects.(subject).scatterArgs = scatterArgs;
            
            % bout size
            bout_size_vals = resTable.bout_size(il_table_values);
            mean_val = nanmean(bout_size_vals);
            data_bout_size.(group).(slcn).Subjects.(subject).mean = mean_val;
            data_bout_size.(group).(slcn).Subjects.(subject).scatterArgs = scatterArgs;
            
            % Area under curve
            for i_window = 1:n_windows
                window = windows{i_window};
                % Lick rate on 
                vals = AUC.(group).(window).subjects.(subject).(slcn).licks.values;
                vals_mean = nanmean(vals);
                vals_SEM = SEM(vals);
                AUC.(group).(window).subjects.(subject).(slcn).licks.mean = vals_mean;
                AUC.(group).(window).subjects.(subject).(slcn).licks.SEM = vals_SEM;
                AUC.(group).(window).subjects.(subject).(slcn).licks.plotArgs = barArgs;
                
                n_corr = size(Correlations.(group).(window).values_mean, 1);
                Correlations.(group).(window).values_mean(n_corr + 1, 1) = vals_mean;

                % Photometry Zscore
                vals = AUC.(group).(window).subjects.(subject).(slcn).Zscore.values;
                vals_mean = nanmean(vals);
                vals_SEM = SEM(vals);
                AUC.(group).(window).subjects.(subject).(slcn).Zscore.mean = vals_mean;
                AUC.(group).(window).subjects.(subject).(slcn).Zscore.SEM = vals_SEM;
                AUC.(group).(window).subjects.(subject).(slcn).Zscore.plotArgs = barArgs;
                Correlations.(group).(window).values_mean(n_corr + 1, 2) = vals_mean;

                % Photometry Peak   
                vals = AUC.(group).(window).subjects.(subject).(slcn).Z_max.values;
                vals_mean = nanmean(vals);
                vals_SEM = SEM(vals);
                AUC.(group).(window).subjects.(subject).(slcn).Z_max.mean = vals_mean;
                AUC.(group).(window).subjects.(subject).(slcn).Z_max.SEM = vals_SEM;
                AUC.(group).(window).subjects.(subject).(slcn).Z_max.plotArgs = barArgs;
                Correlations.(group).(window).values_mean(n_corr + 1, 3) = vals_mean;

                % Photometry Latency to Peak
                vals = AUC.(group).(window).subjects.(subject).(slcn).Z_peak_lat.values;
                vals_mean = nanmean(vals);
                vals_SEM = SEM(vals);
                AUC.(group).(window).subjects.(subject).(slcn).Z_peak_lat.mean = vals_mean;
                AUC.(group).(window).subjects.(subject).(slcn).Z_peak_lat.SEM = vals_SEM;
                AUC.(group).(window).subjects.(subject).(slcn).Z_peak_lat.plotArgs = barArgs;
                Correlations.(group).(window).values_mean(n_corr + 1, 4) = vals_mean;
                
             end
            
        end
        % photo
        % dF
        vals_dF = resTable.dF(il_table_group_slcn);
        all_vals = cat(2, vals_dF{:});
        vals_mean = nanmean(all_vals, 2);
        vals_SEM = nanstd(all_vals, 0, 2) / sqrt(size(all_vals, 1));
        data_dF.(group).(slcn).mean = vals_mean;
        data_dF.(group).(slcn).SEM = vals_SEM;
        data_dF.(group).(slcn).plotArgs = plotArgs;
        
        % Zscore
        vals_Zscore = resTable.Zscore(il_table_group_slcn);
        all_vals = cat(2, vals_Zscore{:});
        vals_mean = nanmean(all_vals, 2);
        vals_SEM = nanstd(all_vals, 0, 2) / sqrt(size(all_vals, 1));
        data_Zscore.(group).(slcn).mean = vals_mean;
        data_Zscore.(group).(slcn).SEM = vals_SEM;
        data_Zscore.(group).(slcn).plotArgs = plotArgs;
        
        %WINDOWS
            %Z-score
            for i_window = 1:n_windows
                window = windows{i_window};
                
                vals = AUC.(group).(window).(slcn).Zscore.values;
                vals_mean = nanmean(vals);
                vals_SEM = SEM(vals);
                AUC.(group).(window).(slcn).Zscore.mean = vals_mean;
                AUC.(group).(window).(slcn).Zscore.SEM = vals_SEM;
                AUC.(group).(window).(slcn).Zscore.plotArgs = barArgs;
                %Area Under Curve
                vals = AUC.(group).(window).(slcn).licks.values;
                vals_mean = nanmean(vals);
                vals_SEM = SEM(vals);
                AUC.(group).(window).(slcn).licks.mean = vals_mean;
                AUC.(group).(window).(slcn).licks.SEM = vals_SEM;
                AUC.(group).(window).(slcn).licks.plotArgs = barArgs;
                %Peak (maximum)
                vals = AUC.(group).(window).(slcn).Z_max.values;
                vals_mean = nanmean(vals);
                vals_SEM = SEM(vals);
                AUC.(group).(window).(slcn).Z_max.mean = vals_mean;
                AUC.(group).(window).(slcn).Z_max.SEM = vals_SEM;
                AUC.(group).(window).(slcn).Z_max.plotArgs = barArgs;
                %Peak latency
                vals = AUC.(group).(window).(slcn).Z_peak_lat.values;
                vals_mean = nanmean(vals);
                vals_SEM = SEM(vals);
                AUC.(group).(window).(slcn).Z_peak_lat.mean = vals_mean;
                AUC.(group).(window).(slcn).Z_peak_lat.SEM = vals_SEM;
                AUC.(group).(window).(slcn).Z_peak_lat.plotArgs = barArgs;

            end
    end    
end

% Pearson correlationships for mean values
for i_group = 1:n_groups
    group = groups{i_group};
    
    for i_window = 1:n_windows
        window = windows{i_window};
        [rho, pval] = corr(Correlations.(group).(window).values_mean);
        Correlations.(group).(window).Results_mean.AUC_corr = rho(1,2);
        Correlations.(group).(window).Results_mean.AUC_p = pval(1,2);
        Correlations.(group).(window).Text_mean.r_AUC = strcat('r = ', num2str(rho(1,2)), ', p =', num2str(pval(1,2)));
        
        Correlations.(group).(window).Results_mean.peak_corr = rho(1,3);
        Correlations.(group).(window).Results_mean.peak_p = pval(1,3);
        Correlations.(group).(window).Text_mean.r_peak = strcat('r = ', num2str(rho(1,3)), ', p =', num2str(pval(1,3)));
        
        Correlations.(group).(window).Results_mean.latency_corr = rho(1,4);
        Correlations.(group).(window).Results_mean.latency_p = pval(1,4);
        Correlations.(group).(window).Text_mean.r_latency = strcat('r = ', num2str(rho(1,4)), ', p =', num2str(pval(1,4)));
           
    end % windows
    
end % groups 

% legends for mean correlations
for i_window = 1:n_windows
    window = windows{i_window};
    
    for i_group = 1:n_groups
        group = groups{i_group};
        
        Correlations.Text.(window).Text_AUC_mean{i_group} =  Correlations.(group).(window).Text_mean.r_AUC;
        Correlations.Text.(window).Text_peak_mean{i_group} =  Correlations.(group).(window).Text_mean.r_peak; 
        Correlations.Text.(window).Text_latency_mean{i_group} =  Correlations.(group).(window).Text_mean.r_latency;
        
        % Fit Line for Pearson correlation (Mean)
        % get data mean points
        x = Correlations.(group).(window).values_mean(:,1);
        y_Zs_mean = Correlations.(group).(window).values_mean(:,2);
        y_peak_mean = Correlations.(group).(window).values_mean(:,3);
        y_lat_mean = Correlations.(group).(window).values_mean(:,4);
        % Save beta coefficients for line
        Correlations.(group).(window).fitline.Zs_mean = polyfit(x, y_Zs_mean, 1);
        Correlations.(group).(window).fitline.peak_mean = polyfit(x, y_peak_mean, 1);
        Correlations.(group).(window).fitline.lat_mean = polyfit(x, y_lat_mean, 1);
        
    end
    
end

%% Plotting

% plot
for i_group = 1:n_groups
    group = groups{i_group};
    il_table_group = strcmp(resTable.Group, group);
    subjects_group = unique(resTable.Subject(il_table_group));
    
    % initialize figures
    fig_name = strcat(group, ' - lick rate');
    f_lick_rate = figure('NumberTitle', 'off', 'Name', fig_name);
    hold on
    fig_name = strcat(group, ' - bout size');
    f_bout_size = figure('NumberTitle', 'off', 'Name', fig_name);
    hold on
    
    % iterate over solutions
    for i_slcn = 1:n_slcns
        slcn = slcn_seq{i_slcn};
        x = i_slcn;
        
        % plot lick rate
        y = data_lick_rate.(group).(slcn).mean;
        sem_error = data_lick_rate.(group).(slcn).SEM;
        barArgs = data_lick_rate.(group).(slcn).barArgs;
        
        figure(f_lick_rate);
        bar(x, y, barArgs);
        errorbar(x, y, sem_error, 'Color', 'k');
        title(' Mean Lick Rate for all sessions');
        xlabel('Sucrose %');
        ylabel('Licks/s');
        set(gca, 'xtick', 1:length(bar_x_label), 'xticklabel', (bar_x_label));
        
        % plot bout size
        y = data_bout_size.(group).(slcn).mean;
        sem_error = data_lick_rate.(group).(slcn).SEM;
        barArgs = data_bout_size.(group).(slcn).barArgs;
        
        figure(f_bout_size);
        bar(x, y, barArgs);
        errorbar(x, y, sem_error, 'Color', 'k');
        title(' Mean Bout Size for all sessions');
        xlabel('Sucrose %');
        ylabel('Bout size (s)');
        set(gca, 'xtick', 1:length(bar_x_label), 'xticklabel',(bar_x_label));
        
        % plot markers for subjects in lick rate and bout size
        for i_subject = 1:numel(subjects_group)
            subject = subjects_group{i_subject};
            
            % plot lick rate subject markers
            y =  data_lick_rate.(group).(slcn).Subjects.(subject).mean;
            scatterArgs = data_lick_rate.(group).(slcn).Subjects.(subject).scatterArgs;
            figure(f_lick_rate);
            h = scatter(x, y, marker_size);
            set(h, scatterArgs);
            
            % plot bout size subject markers
            y =  data_bout_size.(group).(slcn).Subjects.(subject).mean;
            scatterArgs = data_bout_size.(group).(slcn).Subjects.(subject).scatterArgs;
            figure(f_bout_size);
            h = scatter(x, y, marker_size);
            set(h, scatterArgs);           
            
        end
    end
end

fig_name = strcat(group, ' - Photometry');
f_photo = figure('NumberTitle', 'off', 'Name', fig_name);

axlim_licks.xlim = [lags(1) lags(end)];
axlim_licks.ylim = photo_licks_y_lim;

axlim_Zscore.xlim = [photo_x_time(1) photo_x_time(end)];
axlim_Zscore.ylim = photo_Zscore_y_lim;
x_patch = [0, 0, 4, 4];
y_patch = [-5, 15, 15, -5];


for i_slcn = 1:n_slcns
    slcn = slcn_seq{i_slcn};
    shadeArgs.FaceColor = color_slcn.(slcn);
    shadeArgs.EdgeColor = 'none';
    shadeArgs.FaceAlpha = [0.3];
    
    for i_group = 1:n_groups
        group = groups{i_group};
        % plot Z-score
        % ax properties
        subplot(2, n_slcns, i_slcn);
        set(gca, axlim_Zscore);
        xlabel('Time');
        ylabel('Z-score');
        hold on
        % get data
        x = photo_x_time;
        Zscore_y = data_Zscore.(group).(slcn).mean;
        y = smoothPSTH(Zscore_y, args.L);
        yy_line = get(gca, 'ylim');
        sem_error = data_Zscore.(group).(slcn).SEM;
        plotArgs.Color = color_group.(group);
        plotArgs.LineStyle = '-';
        % plot
        s = patch(x_patch, y_patch, 'w');
        set(s, shadeArgs);
        plot(x, y, plotArgs);
        shadedErrorBar(x, y, sem_error, plotArgs, 1);

        
        % plot lick rate
        % ax properties
        subplot(2, n_slcns, i_slcn + n_slcns);
        set(gca, axlim_licks);
        xlabel('Time');
        ylabel('Licks/s');
        hold on
        % data
        x = lags;
        y = data_lick_rate.(group).(slcn).PSTH.mean;
        yy_line = get(gca, 'ylim');
        sem_error = data_lick_rate.(group).(slcn).PSTH.SEM;
        plotArgs.Color = color_group.(group);
        plotArgs.LineStyle = '-';

        % plot
        s = patch(x_patch, y_patch, 'w');
        set(s, shadeArgs);
        plot(x, y, plotArgs);
        shadedErrorBar(x, y, sem_error, plotArgs, 0);
        
    end
end
ax1 = subplot(2, n_slcns, 1);
legend(ax1, groups');

ax2 = subplot(2, n_slcns, n_slcns+1);
legend(ax2, groups')

%% plot Z-score & PSTH licks for every solution

for i_group = 1:n_groups
    group = groups{i_group};
    
    fig_name = strcat(group, ' - Zscore-Photometry');
    figure('NumberTitle', 'off', 'Name', fig_name);

    for i_slcn = 1:n_slcns
        slcn = slcn_seq{i_slcn};
        
        %plot photometry Zscore
        subplot(2, 1, 1);
        hold on
        xlabel('Time');
        ylabel('Z-score');
       
        % get data
        x = photo_x_time;
        Zscore_y = data_Zscore.(group).(slcn).mean;
        y = smoothPSTH(Zscore_y, args.L);
        sem_error = data_Zscore.(group).(slcn).SEM;
        plotArgs = data_Zscore.(group).(slcn).plotArgs;
        % plot
        plot(x, y, plotArgs);
        shadedErrorBar(x, y, sem_error, plotArgs, 1);
        
        % plot lick rate
        % ax properties
        subplot(2, 1, 2);
        hold on
        set(gca, axlim_licks);
        xlabel('Time');
        ylabel('Licks/s');
        
        % data
        x = lags;
        y = data_lick_rate.(group).(slcn).PSTH.mean;
        sem_error = data_lick_rate.(group).(slcn).PSTH.SEM;
        plotArgs = data_lick_rate.(group).(slcn).PSTH.plotArgs;
        % plot
        plot(x, y, plotArgs);
        shadedErrorBar(x, y, sem_error, plotArgs, 1);
        
    end
end


figure('NumberTitle', 'off', 'Name', 'Photometry Z-score difference');
subplot(1, 1, 1);
hold on
xlabel('Time (s)');
ylabel('Z-score');
x = photo_x_time;

il_non_ref_groups = cellfun(@(g) (~strcmp(g, reference_group)), groups);
non_ref_groups = groups(il_non_ref_groups);
for i_group = 1:length(non_ref_groups)
    group = non_ref_groups{i_group};
    
    for i_slcn = 1:n_slcns
        slcn = slcn_seq{i_slcn};
        %
        y_group = data_Zscore.(group).(slcn).mean;
        y_ref = data_Zscore.(reference_group).(slcn).mean;
        y = y_group - y_ref;
        y = smoothPSTH(y, args.L);
%         sem_error = data_lick_rate.(group).(slcn).PSTH.SEM;
        plotArgs = data_lick_rate.(group).(slcn).PSTH.plotArgs;

        % plot
        plot(x, y, plotArgs);
%         shadedErrorBar(x, y, sem_error, plotArgs, 1);
    end
end


%% MODIFIQUE AQUI -Area under curve for Zscore & lick rate-

for i_window = 1:n_windows
    window = windows{i_window};


    for i_group = 1:n_groups
        group = groups{i_group};
        subjects = fieldnames(AUC.(group).(window).subjects);
        
        fig_name = strcat('AUC & Correlations Window -', num2str(i_window));
        figure('NumberTitle', 'off', 'Name', fig_name);
        
        for i_slcn = 1:n_slcns
            slcn = slcn_seq{i_slcn};
            x = i_slcn;


            if i_group == 2;
               scatterArgs.markerFaceColor = 'none';
               scatterArgs.markerEdgeColor = color_slcn.(slcn);

            else
               scatterArgs.markerFaceColor = color_slcn.(slcn);
               scatterArgs.markerEdgeColor = color_slcn.(slcn);
            end

            % positive area under curve for lick rate PSTH
            y = AUC.(group).(window).(slcn).licks.mean;
            sem_error = AUC.(group).(window).(slcn).licks.SEM;
            barArgs = AUC.(group).(window).(slcn).licks.plotArgs;
            % plot
            subplot(5, 4, 1);
            hold on
            bar(x, y, barArgs);
            errorbar(x, y, sem_error, 'Color', 'k');
            xlabel('Sucrose %');
            ylabel('AUC lick rate');
            set(gca, 'xtick', 1:length(bar_x_label), 'xticklabel',(bar_x_label));
            for i_s = 1:numel(subjects);
                subject = subjects{i_s};
                scatterArgs.Marker = marker_seq{i_s};
                y = AUC.(group).(window).subjects.(subject).(slcn).licks.mean;
                h = scatter(x, y, marker_size);
                set(h, scatterArgs);
            end

            % positive area under curve for Z-score
            y = AUC.(group).(window).(slcn).Zscore.mean;
            sem_error = AUC.(group).(window).(slcn).Zscore.SEM;
            barArgs = AUC.(group).(window).(slcn).Zscore.plotArgs;
            % plot
            subplot(5, 4, 2);
            hold on
            bar(x, y, barArgs);
            errorbar(x, y, sem_error, 'Color', 'k');
            xlabel('Sucrose %');
            ylabel('AUC Z-score');
            set(gca, 'xtick', 1:length(bar_x_label), 'xticklabel',(bar_x_label));
            for i_s = 1:numel(subjects);
                subject = subjects{i_s};
                scatterArgs.Marker = marker_seq{i_s};
                y = AUC.(group).(window).subjects.(subject).(slcn).Zscore.mean;
                h = scatter(x, y, marker_size);
                set(h, scatterArgs);
            end

            % Z-score maximum
            y = AUC.(group).(window).(slcn).Z_max.mean;
            sem_error = AUC.(group).(window).(slcn).Z_max.SEM;
            barArgs = AUC.(group).(window).(slcn).Z_max.plotArgs;
            % plot
            subplot(5, 4, 3);
            hold on
            bar(x, y, barArgs);
            errorbar(x, y, sem_error, 'Color', 'k');
            xlabel('Sucrose %');
            ylabel('Z-score max');
            set(gca, 'xtick', 1:length(bar_x_label), 'xticklabel',(bar_x_label));
            for i_s = 1:numel(subjects);
                subject = subjects{i_s};
                scatterArgs.Marker = marker_seq{i_s};
                y = AUC.(group).(window).subjects.(subject).(slcn).Z_max.mean;
                h = scatter(x, y, marker_size);
                set(h, scatterArgs);
            end

            % Z-score peak latency
            y = AUC.(group).(window).(slcn).Z_peak_lat.mean;
            sem_error = AUC.(group).(window).(slcn).Z_peak_lat.SEM;
            barArgs = AUC.(group).(window).(slcn).Z_peak_lat.plotArgs;
            % plot
            subplot(5, 4, 4);
            hold on
            bar(x, y, barArgs);
            errorbar(x, y, sem_error, 'Color', 'k');
            xlabel('Sucrose %');
            ylabel('Z-score peak latency');
            set(gca, 'xtick', 1:length(bar_x_label), 'xticklabel',(bar_x_label));
            for i_s = 1:numel(subjects);
                subject = subjects{i_s};
                scatterArgs.Marker = marker_seq{i_s};
                y = AUC.(group).(window).subjects.(subject).(slcn).Z_peak_lat.mean;
                h = scatter(x, y, marker_size);
                set(h, scatterArgs);
            end

            % relationship of positive areas under curve       
            ax1 = subplot(5, 3, [4 7]);
            hold on
            xlabel('AUC Lick rate');
            ylabel('AUC Zscore');
            
            ax2 = subplot(5, 3, [5 8]);
            hold on
            xlabel('AUC Lick rate');
            ylabel('Z-score maximum');

            ax3 = subplot(5, 3, [6 9]);
            hold on
            xlabel('AUC Lick rate');
            ylabel('Z-score peak latency');


            for i_s = 1:numel(subjects);
                subject = subjects{i_s};
                scatterArgs.Marker = marker_seq{i_s};
                x = AUC.(group).(window).subjects.(subject).(slcn).licks.values;

                % licks vs Zscore
                y = AUC.(group).(window).subjects.(subject).(slcn).Zscore.values;
                subplot(ax1);
                h = scatter(x, y, marker_size);
                set(h, scatterArgs)              

                % licks vs Z_max
                y = AUC.(group).(window).subjects.(subject).(slcn).Z_max.values;
                subplot(ax2);
                h = scatter(x, y, marker_size);
                set(h, scatterArgs)
                

                % licks vs Z_peak_lat
                y = AUC.(group).(window).subjects.(subject).(slcn).Z_peak_lat.values;
                subplot(ax3);
                h = scatter(x, y, marker_size);
                set(h, scatterArgs)
                

            end
            
            
            % Correlation per subject whit means
            ax4 = subplot(5, 3, [10 13]);
            hold on
            xlabel('AUC Lick rate');
            ylabel('AUC Zscore');

            ax5 = subplot(5, 3, [11 14]);
            hold on
            xlabel('AUC Lick rate');
            ylabel('Z-score maximum');

            ax6 = subplot(5, 3, [12 15]);
            hold on
            xlabel('AUC Lick rate');
            ylabel('Z-score peak latency');


            for i_s = 1:numel(subjects);
                subject = subjects{i_s};
                scatterArgs.Marker = marker_seq{i_s};
                x = AUC.(group).(window).subjects.(subject).(slcn).licks.mean;
                

                % licks vs Zscore
                y = AUC.(group).(window).subjects.(subject).(slcn).Zscore.mean;
                subplot(ax4);
                h = scatter(x, y, marker_size);
                set(h, scatterArgs)
                

                % licks vs Z_max
                y = AUC.(group).(window).subjects.(subject).(slcn).Z_max.mean;
                subplot(ax5);
                h = scatter(x, y, marker_size);
                set(h, scatterArgs)
                

                % licks vs Z_peak_lat
                y = AUC.(group).(window).subjects.(subject).(slcn).Z_peak_lat.mean;
                subplot(ax6);
                h = scatter(x, y, marker_size);
                set(h, scatterArgs)
                

            end
            
        end %solutions
        
        
            % Plot fit line for all data      
            ln = refline(ax1, Correlations.(group).(window).fitline.Zs);
            set(ln, 'Color', 'k');
            ln = refline(ax2, Correlations.(group).(window).fitline.peak);
            set(ln, 'Color', 'k');
            ln = refline(ax3, Correlations.(group).(window).fitline.lat);
            set(ln, 'Color', 'k');
            % Plot fit line for mean data
            ln = refline(ax4, Correlations.(group).(window).fitline.Zs_mean);
            set(ln, 'Color', 'k');
            ln = refline(ax5, Correlations.(group).(window).fitline.peak_mean);
            set(ln, 'Color', 'k');
            ln = refline(ax6, Correlations.(group).(window).fitline.lat_mean);
            set(ln, 'Color', 'k');
        
    end %groups
    
    legend(ax1, Correlations.Text.(window).Text_AUC);

    
    legend(ax2, Correlations.Text.(window).Text_peak);
    legend(ax3, Correlations.Text.(window).Text_latency);
    legend(ax4, Correlations.Text.(window).Text_AUC_mean);
    legend(ax5, Correlations.Text.(window).Text_peak_mean);
    legend(ax6, Correlations.Text.(window).Text_latency_mean);
    
    
end %windows



%% Area under curve for Zscore & lick rate 
%     
% for i_group = 1:n_groups
%     group = groups{i_group};
%     fig_name = strcat('AUC - ', group);
%     figure('NumberTitle', 'off', 'Name', fig_name);
% 
%     subjects = fieldnames(AUC.(group).W1.subjects);
%     
%     for i_slcn = 1:n_slcns
%         slcn = slcn_seq{i_slcn};
%         x = i_slcn;
% 
%         
%         if i_group == 2;
%            scatterArgs.markerFaceColor = 'none';
%            scatterArgs.markerEdgeColor = color_slcn.(slcn);
%         
%         else
%            scatterArgs.markerFaceColor = color_slcn.(slcn);
%            scatterArgs.markerEdgeColor = color_slcn.(slcn);
%         end
%         
%         % positive area under curve for lick rate PSTH
%         y = AUC.(group).(slcn).licks.mean;
%         sem_error = AUC.(group).(slcn).licks.SEM;
%         barArgs = AUC.(group).(slcn).licks.plotArgs;
%         % plot
%         subplot(5, 4, 1);
%         hold on
%         bar(x, y, barArgs);
%         errorbar(x, y, sem_error, 'Color', 'k');
%         xlabel('Sucrose %');
%         ylabel('AUC lick rate');
%         set(gca, 'xtick', 1:length(bar_x_label), 'xticklabel',(bar_x_label));
%         for i_s = 1:numel(subjects);
%             subject = subjects{i_s};
%             scatterArgs.Marker = marker_seq{i_s};
%             y = AUC.(group).subjects.(subject).(slcn).licks.mean;
%             h = scatter(x, y, marker_size);
%             set(h, scatterArgs);
%         end
%         
%         % positive area under curve for Z-score
%         y = AUC.(group).(slcn).Zscore.mean;
%         sem_error = AUC.(group).(slcn).Zscore.SEM;
%         barArgs = AUC.(group).(slcn).Zscore.plotArgs;
%         % plot
%         subplot(5, 4, 2);
%         hold on
%         bar(x, y, barArgs);
%         errorbar(x, y, sem_error, 'Color', 'k');
%         xlabel('Sucrose %');
%         ylabel('AUC Z-score');
%         set(gca, 'xtick', 1:length(bar_x_label), 'xticklabel',(bar_x_label));
%         for i_s = 1:numel(subjects);
%             subject = subjects{i_s};
%             scatterArgs.Marker = marker_seq{i_s};
%             y = AUC.(group).subjects.(subject).(slcn).Zscore.mean;
%             h = scatter(x, y, marker_size);
%             set(h, scatterArgs);
%         end
%         
%         % Z-score maximum
%         y = AUC.(group).(slcn).Z_max.mean;
%         sem_error = AUC.(group).(slcn).Z_max.SEM;
%         barArgs = AUC.(group).(slcn).Z_max.plotArgs;
%         % plot
%         subplot(5, 4, 3);
%         hold on
%         bar(x, y, barArgs);
%         errorbar(x, y, sem_error, 'Color', 'k');
%         xlabel('Sucrose %');
%         ylabel('Z-score max');
%         set(gca, 'xtick', 1:length(bar_x_label), 'xticklabel',(bar_x_label));
%         for i_s = 1:numel(subjects);
%             subject = subjects{i_s};
%             scatterArgs.Marker = marker_seq{i_s};
%             y = AUC.(group).subjects.(subject).(slcn).Z_max.mean;
%             h = scatter(x, y, marker_size);
%             set(h, scatterArgs);
%         end
%         
%         % Z-score peak latency
%         y = AUC.(group).(slcn).Z_peak_lat.mean;
%         sem_error = AUC.(group).(slcn).Z_peak_lat.SEM;
%         barArgs = AUC.(group).(slcn).Z_peak_lat.plotArgs;
%         % plot
%         subplot(5, 4, 4);
%         hold on
%         bar(x, y, barArgs);
%         errorbar(x, y, sem_error, 'Color', 'k');
%         xlabel('Sucrose %');
%         ylabel('Z-score peak latency');
%         set(gca, 'xtick', 1:length(bar_x_label), 'xticklabel',(bar_x_label));
%         for i_s = 1:numel(subjects);
%             subject = subjects{i_s};
%             scatterArgs.Marker = marker_seq{i_s};
%             y = AUC.(group).subjects.(subject).(slcn).Z_peak_lat.mean;
%             h = scatter(x, y, marker_size);
%             set(h, scatterArgs);
%         end
%         
%         % relationship of positive areas under curve
%         
%         ax1 = subplot(5, 3, [4 7]);
%         hold on
% %         title({'AUC Zscore vs', 'AUC Lick rate per Session'});
%         xlabel('AUC Lick rate');
%         ylabel('AUC Zscore');
% 
%         ax2 = subplot(5, 3, [5 8]);
%         hold on
% %         title({'Mean AUC Zscore vs', 'AUC Lick rate'});
%         xlabel('AUC Lick rate');
%         ylabel('Z-score maximum');
% 
%         ax3 = subplot(5, 3, [6 9]);
%         hold on
% %         title({'Mean AUC Zscore vs', 'AUC Lick rate'});
%         xlabel('AUC Lick rate');
%         ylabel('Z-score peak latency');
% 
%         
%         for i_s = 1:numel(subjects);
%             subject = subjects{i_s};
%             scatterArgs.Marker = marker_seq{i_s};
%             x = AUC.(group).subjects.(subject).(slcn).licks.mean;
%             
%             % licks vs Zscore
%             y = AUC.(group).subjects.(subject).(slcn).Zscore.mean;
%             subplot(ax1);
%             h = scatter(x, y, marker_size);
%             set(h, scatterArgs)
%             
%             % licks vs Z_max
%             y = AUC.(group).subjects.(subject).(slcn).Z_max.mean;
%             subplot(ax2);
%             h = scatter(x, y, marker_size);
%             set(h, scatterArgs)
%             
%             % licks vs Z_peak_lat
%             y = AUC.(group).subjects.(subject).(slcn).Z_peak_lat.mean;
%             subplot(ax3);
%             h = scatter(x, y, marker_size);
%             set(h, scatterArgs)
% 
%         end
% %         x = AUC.(group).(slcn).licks.values;
% %         y = AUC.(group).(slcn).Zscore.values;
%         scatterArgs.Marker = marker_seq{i_slcn};
% %         % plot
% %         subplot(5, 3, [4 7]);
% %         hold on
% %         h = scatter(x, y, marker_size);
% %         set(h, scatterArgs)
% %         title({'AUC Zscore vs', 'AUC Lick rate per Session'});
% %         xlabel('AUC Lick rate');
% %         ylabel('AUC Zscore');
% 
%         subplot(5, 3, [10 13]);
%         hold on
%         mean_x = AUC.(group).(slcn).licks.mean;
%         mean_y = AUC.(group).(slcn).Zscore.mean;
%         h = scatter(mean_x, mean_y, marker_size);
%         set(h, scatterArgs)
%         title({'Mean AUC Zscore vs', 'AUC Lick rate'});
%         xlabel('AUC Lick rate');
%         ylabel('AUC Zscore');
%         
%         
%         % relationship of AUC licks vs. Z-score maximum
% %         x = AUC.(group).(slcn).licks.values;
% %         y = AUC.(group).(slcn).Z_max.values;
% %         % plot
% %         subplot(5, 3, [5 8]);
% %         hold on
% %         h = scatter(x, y, marker_size);
% %         set(h, scatterArgs)
% %         title({'AUC Zscore vs', 'AUC Lick rate per Session'});
% %         xlabel('AUC Lick rate');
% %         ylabel('AUC Zscore');
% 
%         subplot(5, 3, [11 14]);
%         hold on
%         mean_x = AUC.(group).(slcn).licks.mean;
%         mean_y = AUC.(group).(slcn).Z_max.mean;
%         h = scatter(mean_x, mean_y, marker_size);
%         set(h, scatterArgs)
%         title({'Mean AUC Zscore vs', 'Z-score maximum'});
%         xlabel('AUC Lick rate');
%         ylabel('Z-score maximum');
%         
%         % relationship of positive areas under curve
% %         x = AUC.(group).(slcn).licks.values;
% %         y = AUC.(group).(slcn).Z_peak_lat.values;
% %         % plot
% %         subplot(5, 3, [6 9]);
% %         hold on
% %         h = scatter(x, y, marker_size);
% %         set(h, scatterArgs)
% %         title({'AUC Zscore vs', 'AUC Lick rate per Session'});
% %         xlabel('AUC Lick rate');
% %         ylabel('AUC Zscore');
% 
%         subplot(5, 3, [12 15]);
%         hold on
%         mean_x = AUC.(group).(slcn).licks.mean;
%         mean_y = AUC.(group).(slcn).Z_peak_lat.mean;
%         h = scatter(mean_x, mean_y, marker_size);
%         set(h, scatterArgs)
%         title({'Mean AUC Zscore vs', 'Z-score peak latency'});
%         xlabel('AUC Lick rate');
%         ylabel('Z-score peak latency');
%         
%         
%     end
% end

