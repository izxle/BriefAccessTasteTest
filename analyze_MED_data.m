function data = analyze_MED_data(varargin)
%% AUTHOR    : Oscar X. Guerrero-Gutierrez
%% $DATE     : 01-Feb-2019 $
%% DEVELOPED : (R2015a)
%% FILENAME  : analyze_MED_data.m
%% Parameters
% argument parser
pArgs = inputParser;
pMED = inputParser;
% parameter definition
attrMED = {'positive', 'increasing', 'nonempty'};
validateMEDvalues = @(x) validateattributes(x, {'double'}, attrMED);
validateValves = @(x) cellfun(validateMEDvalues, struct2cell(x));
% required parameters
% TODO: use validate functions
pArgs.addParameter('data_MED', 0);
pArgs.addParameter('name_table', 0);
% aliases for key values
pArgs.addParameter('events_label', 'Events');
pArgs.addParameter('licks_label', 'Licks');
% analysis parameters
pArgs.addParameter('rwd_period', 4);
pArgs.addParameter('max_IBI', 0.5);  % Inter-Bout Interval threshold
% window around event of trial
pArgs.addParameter('minlag', -3);
pArgs.addParameter('maxlag', 6);
% pSTH parameters
pArgs.addParameter('resolution', 0.05);
pArgs.addParameter('L', 3);  % smooth factor

% parse arguments
pArgs.parse(varargin{:});
args = pArgs.Results;

% pMED.addParameter('events', 0, validateMEDvalues);
% pMED.addParameter('licks', 0, validateMEDvalues);
% pMED.addParameter('valves', 0, validateValves);
% pMED.addParameter('rwd', 0, validateMEDvalues);
% %
% pMED.parse(args.values_MED);
% values_MED = pMED.Results;
% TODO: check name_table and values_MED
%% Parse raw MED data values
data = parse_MED_values(args.data_MED, args.name_table);

%% initialize variables
solutions = data.Solutions;
events_all = data.Values.(args.events_label);
licks = data.Values.(args.licks_label);
interlick_interval = diff(licks);
il_licks_bout_init = [1; interlick_interval >= args.max_IBI];
i_licks_bout_init = find(il_licks_bout_init);
n_events = numel(events_all);
events_slcn_type = zeros(n_events, 1);
lags = lagsfun(args.minlag, args.maxlag, args.resolution);


slcn_names = fields(solutions);
cell_init = cell(n_events, 1);
nan_init = nan(n_events, 1);

% initialize trials table
trialsTable = table;
trialsTable.solution = cell_init;
trialsTable.values = cell_init;
trialsTable.type = cell_init;
trialsTable.event = nan_init;
trialsTable.first_bout = cell_init;
trialsTable.bout_size = nan_init;
trialsTable.licks = cell_init;
trialsTable.lick_rate = nan_init;
trialsTable.licks_out = cell_init;
trialsTable.licks_after = cell_init;
trialsTable.licks_after_len = nan_init;

%% main

% iterate over solutions
for i_slcn = 1:numel(slcn_names)
    slcn = slcn_names{i_slcn};
    % round to bypass float uncertainty
    values = round(solutions.(slcn).values, 2);
    i_events_slcn = find(ismember(events_all, values));
    events_slcn_type(i_events_slcn) = i_slcn;
    
    events_slcn = events_all(i_events_slcn);
    n_trials_slcn = length(events_slcn);
    % initialize slcnTable
    cell_init = cell(n_trials_slcn, 1);
    nan_init = nan(n_trials_slcn, 1);
    slcnTable = table;
    slcnTable.index = nan_init;
    slcnTable.values = cell_init;
    slcnTable.type = cell_init;
    slcnTable.event = nan_init;
    slcnTable.first_bout = cell_init;
    slcnTable.bout_size = nan_init;
    slcnTable.licks = cell_init;
    slcnTable.lick_rate = nan_init;
    slcnTable.licks_out = cell_init;
    slcnTable.licks_after = cell_init;
    slcnTable.licks_after_len = nan_init;
    
    % iterate over trials of one solution
    for i_trial_slcn = 1:n_trials_slcn
        i_trial = i_events_slcn(i_trial_slcn);
        init = events_slcn(i_trial_slcn);
        fin = init + args.rwd_period;
        % get values inside init-fin window
        window_trial = (init <= values) & (values < fin);
        values_trial = round(values(window_trial), 2);
        
        init_window = args.minlag + init;
        fin_window = init + args.maxlag;
        il_licks = (init_window <= licks) & (licks < fin_window);
        licks_trial = round(licks(il_licks), 2);
        
        % get licks out of trial
        il_licks_out = ~ismember(licks_trial, values_trial);
        licks_out = licks_trial(il_licks_out);
        if isempty(licks_out)
            licks_after = nan;
            licks_after_len = nan;
        else
            licks_out = licks_out - init;  % relative values
            % get only licks after trial for length
            licks_after = licks_out(licks_out > 0);
            if isempty(licks_after)
                licks_after_len = nan;
            else
                licks_after_len = licks_after(end) - args.rwd_period;
            end
        end
        
        % get first bout of licks
        i_first_lick = find(licks == init);
        try
        i_bout = find(i_licks_bout_init == i_first_lick);
        catch
            a
        end
        i_licks_bout = i_licks_bout_init(i_bout);
        if i_bout ~= length(i_licks_bout_init)  % check for last bout
            i_licks_second_bout = i_licks_bout_init(i_bout + 1);
            bout = licks(i_licks_bout:i_licks_second_bout-1);
        else
            bout = licks(i_licks_bout:end);
        end
        % check for empty bout
        if isempty(bout) || length(bout) == 1
            bout_size = 0;
            lick_rate = 0;
        else
            bout_size = bout(end) - bout(1);
            lick_rate = length(bout) / args.rwd_period;
        end
        % label trial type
        if bout_size >= args.rwd_period - args.max_IBI
            event_type = 'complete';
        else
            event_type = 'incomplete';
        end
        
        slcnTable.index(i_trial_slcn) = i_trial;
        slcnTable.values{i_trial_slcn} = values_trial;
        slcnTable.first_bout{i_trial_slcn} = bout;
        slcnTable.type{i_trial_slcn} = event_type;
        slcnTable.bout_size(i_trial_slcn) = bout_size;
        slcnTable.lick_rate(i_trial_slcn) = lick_rate;
        slcnTable.event(i_trial_slcn) = init;  % event: first_lick
        slcnTable.licks{i_trial_slcn} = licks_trial;
        slcnTable.licks_out{i_trial_slcn} = licks_out;
        slcnTable.licks_after{i_trial_slcn} = licks_after;
        slcnTable.licks_after_len(i_trial_slcn) = licks_after_len;
        
        
        trialsTable.solution{i_trial} = slcn;
        trialsTable.values{i_trial} = values_trial;
        trialsTable.first_bout{i_trial} = bout;
        trialsTable.type{i_trial} = event_type;
        trialsTable.bout_size(i_trial) = bout_size;
        trialsTable.lick_rate(i_trial) = lick_rate;
        trialsTable.event(i_trial) = init;  % event: first_lick
        trialsTable.licks{i_trial} = licks_trial;
        trialsTable.licks_out{i_trial} = licks_out;
        trialsTable.licks_after{i_trial} = licks_after;
        trialsTable.licks_after_len(i_trial) = licks_after_len;
        
    end  % for trials of slcn
    
    % PSTH
    slcn_licks_all = slcnTable.licks;
    slcn_licks = cat(1, slcn_licks_all{:});
    slcn_events = slcnTable.event;

    [n, ~, ~, ~] = PSTH2(slcn_licks, slcn_events, args.minlag, args.maxlag, args.resolution);
    PSTH = smoothPSTH(n, args.L, lags);
    pPSTH = PSTH;
    pPSTH(PSTH < 0) = 0;
    
    % Window 1
    init = 0;
    fin = 4;
    mask = (init < lags) & (lags < fin);
    x = lags(mask);
    y = pPSTH(mask);
    W1_AUC = trapz(x, y);
    
    % Window 2
    init = 0;
    fin = 2;
    mask = (init < lags) & (lags < fin);
    x = lags(mask);
    y = pPSTH(mask);
    W2_AUC = trapz(x, y);
        
    % Window 3
    init = 2;
    fin = 4;
    mask = (init < lags) & (lags < fin);
    x = lags(mask);
    y = pPSTH(mask);
    W3_AUC = trapz(x, y);
    
        
    % store values
    data.Values.Valves.(slcn).trials = slcnTable; 
    data.Solutions.(slcn).trials = slcnTable;
    data.Solutions.(slcn).PSTH.licks.values = PSTH;
    data.Solutions.(slcn).PSTH.licks.W1.AUC = W1_AUC;
    data.Solutions.(slcn).PSTH.licks.W2.AUC = W2_AUC;
    data.Solutions.(slcn).PSTH.licks.W3.AUC = W3_AUC;
    
end  % for solutions

% store all trial info
data.Trials = trialsTable;
end