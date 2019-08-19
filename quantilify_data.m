function dataTable = quantilify_data(varargin)
% add quantile label to trials table in each session
%% Parameters
pArgs = inputParser;
pArgs.addRequired('dataTable');
%
pArgs.addParameter('p', [0.25 0.75]);
%
pArgs.parse(varargin{:});
args = pArgs.Results;
% aliases
dataTable = args.dataTable;
p = args.p;

%% Main

n_data = size(dataTable, 1);
for i = 1:n_data
    dataStruct = dataTable.Data{i,1};
    trialsTable = dataStruct.Trials;
    
    licks = trialsTable.licks;
    n_licks = cellfun(@length, licks);
    
    qs = quantile(n_licks, p);
    q1 = qs(1); q2 = qs(2);
    
    il_q1 = (n_licks <= q1);
    il_q2 = (q2 <= n_licks);
    
    label = cell(size(licks));
    label(il_q1) = {'q1'};
    label(il_q2) = {'q2'};
    
    % update trials table with quantile information
    trialsTable.quantile = label;
    dataStruct.Trials = trialsTable;
    
    % add labels to every solution
    solutions = fieldnames(dataStruct.Solutions);
    for i_slcn = 1:numel(solutions)
        slcn = solutions{i_slcn};
        ind = dataStruct.Solutions.(slcn).trials.index;
        lbls = label(ind);
        dataStruct.Solutions.(slcn).trials.quantile = lbls;
    end
    
    dataTable.Data{i,1} = dataStruct;
end

end