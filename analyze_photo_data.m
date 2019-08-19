function dataMED = analyze_photo_data(varargin)
%% AUTHOR    : Oscar X. Guerrero-Gutierrez
%% $DATE     : 10-Feb-2019 $
%% DEVELOPED : (R2015a)
%% FILENAME  : analyze_photo_data.m
%% Parameters
% argument parser
pArgs = inputParser;
% parameter definition
% required parameters
pArgs.addRequired('dataMED');
pArgs.addRequired('photoTable');

pArgs.addParameter('photo_label', 'RawF_1')
% analysis parameters
pArgs.addParameter('baseline_start', -3);
pArgs.addParameter('baseline_end', -1);
pArgs.addParameter('rwd_period', 4);
pArgs.addParameter('window_start', -2);
pArgs.addParameter('window_end', 6);
pArgs.addParameter('frame_rate', 30)  % frame rate in Hz
pArgs.addParameter('L', 5);

%WINDOWS
pArgs.addParameter('W2_init', 0);
pArgs.addParameter('W2_fin', 2);
pArgs.addParameter('W3_init', 2);
pArgs.addParameter('W3_fin', 4);

%
validate = @(x) any(validatestring(x, {'max', 'min'}));
pArgs.addParameter('lat_fun', 'max', validate);

% parse arguments
pArgs.parse(varargin{:});
args = pArgs.Results;
FR = args.frame_rate;
photoTable = args.photoTable;
dataMED = args.dataMED;
lat_fun = args.lat_fun;
W2_init = args.W2_init;
W2_fin = args.W2_fin;
W3_init = args.W3_init;
W3_fin = args.W3_fin;

%% main

time = photoTable.time;
n_frames = FR * time;
values = photoTable.(args.photo_label);

trialsTable = dataMED.Trials;
cell_init = cell(size(trialsTable, 1), 1);
nan_init = nan(size(cell_init));
trialsTable.(args.photo_label) = cell_init;
trialsTable.dF = cell_init;
trialsTable.Zscore = cell_init;
trialsTable.W1_AUC = nan_init;
trialsTable.W1_peak = nan_init;
trialsTable.W1_latency = nan_init;
trialsTable.W2_AUC = nan_init;
trialsTable.W2_peak = nan_init;
trialsTable.W2_latency = nan_init;
trialsTable.W3_AUC = nan_init;
trialsTable.W3_peak = nan_init;
trialsTable.W3_latency = nan_init;
trialsTable.baseline_raw = cell_init;
trialsTable.baseline_dF = cell_init;
trialsTable.after_raw = cell_init;
trialsTable.after_dF = cell_init;
trialsTable.after_Zscore = cell_init;
corr_Mat_W1 = [];
corr_Mat_W2 = [];
corr_Mat_W3 = [];


solutions = fieldnames(dataMED.Solutions);
% group = dataMED.Parameters.Group;
for i_solution = 1:numel(solutions)
    slcn = solutions{i_solution};
    slcnTable = dataMED.Solutions.(slcn).trials;
    n_trials = size(slcnTable, 1);
    cell_init = cell(n_trials, 1);
    nan_init = nan(n_trials, 1);
    % initialize table columns
    slcnTable.(args.photo_label) = cell_init;
    slcnTable.dF = cell_init;
    slcnTable.Zscore = cell_init;

    slcnTable.W1_AUC = nan_init;
    slcnTable.W1_peak = nan_init;
    slcnTable.W1_latency = nan_init;
    slcnTable.W2_AUC = nan_init;
    slcnTable.W2_peak = nan_init;
    slcnTable.W2_latency = nan_init;
    slcnTable.W3_AUC = nan_init;
    slcnTable.W3_peak = nan_init;
    slcnTable.W3_latency = nan_init;

    slcnTable.baseline_raw = cell_init;
    slcnTable.baseline_dF = cell_init;
    slcnTable.after_raw = cell_init;
    slcnTable.after_dF = cell_init;
    slcnTable.after_Zscore = cell_init;
    
    for i_trial = 1:n_trials
        init = slcnTable.event(i_trial);
        
        
        init_window = init + args.window_start;
        fin_window = init + args.window_end;
        il_window = (init_window <= time) & (time < fin_window);

        time_window = time(il_window);
        rel_time = time_window - init;
        values_window = values(il_window);
        
        % get baseline
        init_baseline = init + args.baseline_start;
        fin_baseline = init + args.baseline_end;
        il_baseline = (init_baseline <= time) & (time < fin_baseline);
        baseline = values(il_baseline);
        mean_baseline = mean(baseline);
        
        % values on trial (4s)
        fin_trial = init + args.rwd_period;
        il_rwd = (init <= time) & (time < fin_trial);
        time_trial = time(il_rwd);
        rel_time_trial = time_trial - init;
        values_rwd = values(il_rwd);
        
        % values after trial
        fin_trial = init + args.rwd_period;
        il_after = (fin_trial <= time) & (time <= fin_window);
        values_after = values(il_after);
        
        % dF
        dF = (values_window - mean_baseline) / mean_baseline;
        dF_baseline = (baseline - mean_baseline) / mean_baseline;
        dF_rwd = (values_rwd - mean_baseline) / mean_baseline;
        dF_after = (values_after - mean_baseline) / mean_baseline;
        
        % Z-score
        % in case trial started before having a full baseline window
        if init + args.baseline_start < 0
            n_size = size(values_window);
            dF = nan(n_size);
            Zscore = nan(n_size);
            Zscore_rwd = nan(size(values_rwd));
            Zscore_after = nan(size(values_after));
        else
            mean_dF_baseline = mean(dF_baseline);
            std_dF_baseline = std(dF_baseline);
            Zscore = (dF - mean_dF_baseline) / std_dF_baseline;
            Zscore_rwd = (dF_rwd - mean_dF_baseline) / std_dF_baseline;
            Zscore_after = (dF_after - mean_dF_baseline) / std_dF_baseline;
        end
        
        % Z-score area under curve on trial (4s)
        if strcmp(lat_fun, 'max')
            peak_fun = @max;
        else
            peak_fun = @min;
        end
        
        % Window 1 (0-4 s) - reward period
        Zscore_p = smoothPSTH(Zscore_rwd, args.L, rel_time_trial);
        W1_AUC = trapz(rel_time_trial, Zscore_p);
        [W1_peak, i_peak] = peak_fun(Zscore_p);
        W1_latency = rel_time_trial(i_peak);
        
        % Window 2 (0-2 s) *defined by user
        mask = (W2_init < rel_time_trial) & (rel_time_trial < W2_fin);
        W2_x = rel_time_trial(mask);
        W2_y = Zscore_rwd(mask);
        W2_y = smoothPSTH(W2_y, args.L, W2_x);
        W2_AUC = trapz(W2_x, W2_y);
        [W2_peak, i_peak] = peak_fun(W2_y);
        W2_latency = rel_time_trial(i_peak);
        
        % Window 3 (2-4 s) *defined by user
        mask = (W3_init < rel_time_trial) & (rel_time_trial < W3_fin);
        W3_x = rel_time_trial(mask);
        W3_y = Zscore_rwd(mask);
        W3_y = smoothPSTH(W3_y, args.L, W3_x);
        W3_AUC = trapz(W3_x, W3_y);
        [W3_peak, i_peak] = peak_fun(W2_y);
        W3_latency = rel_time_trial(i_peak);
        
        % save to table
        slcnTable.(args.photo_label){i_trial} = values_window;
        slcnTable.dF{i_trial} = dF;
        slcnTable.Zscore{i_trial} = Zscore;
        slcnTable.W1_AUC(i_trial) = W1_AUC;
        slcnTable.W1_peak(i_trial) = W1_peak;
        slcnTable.W1_latency(i_trial) = W1_latency;
        slcnTable.W2_AUC(i_trial) = W2_AUC;
        slcnTable.W2_peak(i_trial) = W2_peak;
        slcnTable.W2_latency(i_trial) = W2_latency;
        slcnTable.W3_AUC(i_trial) = W3_AUC;
        slcnTable.W3_peak(i_trial) = W3_peak;
        slcnTable.W3_latency(i_trial) = W3_latency;
        slcnTable.baseline_raw{i_trial} = baseline;
        slcnTable.baseline_dF{i_trial} = dF_baseline;
        slcnTable.after_raw{i_trial} = values_after;
        slcnTable.after_dF{i_trial} = dF_after;
        slcnTable.after_Zscore{i_trial} = Zscore_after;
       
        % save to general trials table
        index = slcnTable.index(i_trial);
        trialsTable.(args.photo_label){index} = values_window;
        trialsTable.dF{index} = dF;
        trialsTable.Zscore{index} = Zscore;
        trialsTable.W1_AUC(i_trial) = W1_AUC;
        trialsTable.W1_peak(i_trial) = W1_peak;
        trialsTable.W1_latency(i_trial) = W1_latency;
        trialsTable.W2_AUC(i_trial) = W2_AUC;
        trialsTable.W2_peak(i_trial) = W2_peak;
        trialsTable.W2_latency(i_trial) = W2_latency;
        trialsTable.W3_AUC(i_trial) = W3_AUC;
        trialsTable.W3_peak(i_trial) = W3_peak;
        trialsTable.W3_latency(i_trial) = W3_latency;
        trialsTable.baseline_raw{index} = baseline;
        trialsTable.baseline_dF{index} = dF_baseline;
        trialsTable.after_raw{index} = values_after;
        trialsTable.after_dF{index} = dF_after;
        trialsTable.after_Zscore{index} = Zscore_after;

    end  % for trials
    
    %
%     slcn_Zs = slcnTable.Zscore;
%     slcn_Zs = cat(2, slcn_Zs{:});
%     PSTH = sum(slcn_Zs, 2);
    W1_AUC_mean = nanmean(slcnTable.W1_AUC);
    W1_AUC_std = nanstd(slcnTable.W1_AUC);
    W1_peak_max = max(slcnTable.W1_peak);
    W1_peak_mean = nanmean(slcnTable.W1_peak);
    W1_peak_std = nanstd(slcnTable.W1_peak);
    W1_latency_mean = nanmean(slcnTable.W1_latency);
    W1_latency_std = nanstd(slcnTable.W1_latency);
    
    
    W2_AUC_mean = nanmean(slcnTable.W2_AUC);
    W2_AUC_std = nanstd(slcnTable.W2_AUC);
    W2_peak_max = max(slcnTable.W2_peak);
    W2_peak_mean = nanmean(slcnTable.W2_peak);
    W2_peak_std = nanstd(slcnTable.W2_peak);
    W2_latency_mean = nanmean(slcnTable.W2_latency);
    W2_latency_std = nanstd(slcnTable.W2_latency);
    
    W3_AUC_mean = nanmean(slcnTable.W3_AUC);
    W3_AUC_std = nanstd(slcnTable.W3_AUC);
    W3_peak_max = max(slcnTable.W3_peak);
    W3_peak_mean = nanmean(slcnTable.W3_peak);
    W3_peak_std = nanstd(slcnTable.W3_peak);
    W3_latency_mean = nanmean(slcnTable.W3_latency);
    W3_latency_std = nanstd(slcnTable.W3_latency);
    
    % WINDOW_1
    n = size(corr_Mat_W1, 1);
    corr_Mat_W1(n+1,1) = dataMED.Solutions.(slcn).PSTH.licks.W1.AUC;
    corr_Mat_W1(n+1,2) = W1_AUC_mean;
    corr_Mat_W1(n+1,3) = W1_peak_mean;
    corr_Mat_W1(n+1,4) = W1_latency_mean;
    
    % WINDOW_2
    n = size(corr_Mat_W2, 1);
    corr_Mat_W2(n+1,1) = dataMED.Solutions.(slcn).PSTH.licks.W2.AUC;
    corr_Mat_W2(n+1,2) = W2_AUC_mean;
    corr_Mat_W2(n+1,3) = W2_peak_mean;
    corr_Mat_W2(n+1,4) = W2_latency_mean;
    
    %WINDOW_3
    n = size(corr_Mat_W3, 1);
    corr_Mat_W3(n+1,1) = dataMED.Solutions.(slcn).PSTH.licks.W3.AUC;
    corr_Mat_W3(n+1,2) = W1_AUC_mean;
    corr_Mat_W3(n+1,3) = W1_peak_mean;
    corr_Mat_W3(n+1,4) = W1_latency_mean;
    
    % store values
    dataMED.Values.Valves.(slcn).trials = slcnTable;
    dataMED.Solutions.(slcn).trials = slcnTable;
%     dataMED.Solutions.(slcn).PSTH.Zscore.values = PSTH;
    dataMED.Solutions.(slcn).PSTH.Zscore.W1.AUC_mean = W1_AUC_mean;
    dataMED.Solutions.(slcn).PSTH.Zscore.W1.AUC_std = W1_AUC_std;
    dataMED.Solutions.(slcn).PSTH.Zscore.W1.peak_max = W1_peak_max;
    dataMED.Solutions.(slcn).PSTH.Zscore.W1.peak_mean = W1_peak_mean;
    dataMED.Solutions.(slcn).PSTH.Zscore.W1.peak_std = W1_peak_std;
    dataMED.Solutions.(slcn).PSTH.Zscore.W1.latency_mean = W1_latency_mean;
    dataMED.Solutions.(slcn).PSTH.Zscore.W1.latency_std = W1_latency_std;
    dataMED.Solutions.(slcn).PSTH.Zscore.W2.AUC_mean = W2_AUC_mean;
    dataMED.Solutions.(slcn).PSTH.Zscore.W2.AUC_std = W2_AUC_std;
    dataMED.Solutions.(slcn).PSTH.Zscore.W2.peak_max = W2_peak_max;
    dataMED.Solutions.(slcn).PSTH.Zscore.W2.peak_mean = W2_peak_mean;
    dataMED.Solutions.(slcn).PSTH.Zscore.W2.peak_std = W2_peak_std;
    dataMED.Solutions.(slcn).PSTH.Zscore.W2.latency_mean = W2_latency_mean;
    dataMED.Solutions.(slcn).PSTH.Zscore.W2.latency_std = W2_latency_std;
    dataMED.Solutions.(slcn).PSTH.Zscore.W3.AUC_mean = W3_AUC_mean;
    dataMED.Solutions.(slcn).PSTH.Zscore.W3.AUC_std = W3_AUC_std;
    dataMED.Solutions.(slcn).PSTH.Zscore.W3.peak_max = W3_peak_max;
    dataMED.Solutions.(slcn).PSTH.Zscore.W3.peak_mean = W3_peak_mean;
    dataMED.Solutions.(slcn).PSTH.Zscore.W3.peak_std = W3_peak_std;
    dataMED.Solutions.(slcn).PSTH.Zscore.W3.latency_mean = W3_latency_mean;
    dataMED.Solutions.(slcn).PSTH.Zscore.W3.latency_std = W3_latency_std;
    
    
end  % for solutions

% correlations
% WINDOW_1
[rho, pval] = corr(corr_Mat_W1);
AUC_corr = rho(1,2);
AUC_p = pval(1,2);
peak_corr = rho(1,3);
peak_p = pval(1,3);
latency_corr = rho(1,4);
latency_p = pval(1,4);

dataMED.Correlations.W1.AUC_corr = AUC_corr;
dataMED.Correlations.W1.AUC_p = AUC_p;
dataMED.Correlations.W1.peak_corr = peak_corr;
dataMED.Correlations.W1.peak_p = peak_p;
dataMED.Correlations.W1.latency_corr = latency_corr;
dataMED.Correlations.W1.latency_p = latency_p;

%WINDOW_2
[rho, pval] = corr(corr_Mat_W2);
AUC_corr = rho(1,2);
AUC_p = pval(1,2);
peak_corr = rho(1,3);
peak_p = pval(1,3);
latency_corr = rho(1,4);
latency_p = pval(1,4);

dataMED.Correlations.W2.AUC_corr = AUC_corr;
dataMED.Correlations.W2.AUC_p = AUC_p;
dataMED.Correlations.W2.peak_corr = peak_corr;
dataMED.Correlations.W2.peak_p = peak_p;
dataMED.Correlations.W2.latency_corr = latency_corr;
dataMED.Correlations.W2.latency_p = latency_p;

%WINDOW_3
[rho, pval] = corr(corr_Mat_W3);
AUC_corr = rho(1,2);
AUC_p = pval(1,2);
peak_corr = rho(1,3);
peak_p = pval(1,3);
latency_corr = rho(1,4);
latency_p = pval(1,4);

dataMED.Correlations.W3.AUC_corr = AUC_corr;
dataMED.Correlations.W3.AUC_p = AUC_p;
dataMED.Correlations.W3.peak_corr = peak_corr;
dataMED.Correlations.W3.peak_p = peak_p;
dataMED.Correlations.W3.latency_corr = latency_corr;
dataMED.Correlations.W3.latency_p = latency_p;

% save general trials table to dataMED
dataMED.Trials = trialsTable;

end