function mask=BATT_stratifydata(varargin)
%% Parameters
pArgs = inputParser;
pArgs.addRequired('dataTable');
pArgs.addRequired('slcn_seq');

pArgs.parse(varargin{:});
args = pArgs.Results;

%% Alias
dataTable = args.dataTable;
slcn_seq = args.slcn_seq;

%% Initialize variables

groups = unique(dataTable.Group);
n_files = size(dataTable, 1);
n_slcns = length(slcn_seq);

mask = false(n_files, 1);

%% Statistical calculations

for i_file = 1:n_files
    
    sessionData = dataTable.Data{i_file};
    
    trialsTable = sessionData.Trials;
    slcn_tag = trialsTable.solution;
    lick_rate = trialsTable.lick_rate;
    % bout_size = trialsTable.bout_size;
    
    % sort solutions
    [~, i_tags] = ismember(slcn_tag, slcn_seq);
    [~, idx] = sort(i_tags);
    
    lick_rate = lick_rate(idx);
    slcn_tag = slcn_tag(idx);

    %ANOVA ONE WAY CALCULATIONS
    [p, tbl, stats] = anova1(lick_rate, slcn_tag, 'off');
    if p >= 0.05
        continue
    end
    
    % TURKEY Post Hoc
    compare = multcompare(stats,'CType', 'lsd', 'Display', 'off');
    
    p_water_sac5 = compare(n_slcns-1, 6);
    
    if p_water_sac5 < 0.05;
        mask(i_file) = 1;
    end

end


end