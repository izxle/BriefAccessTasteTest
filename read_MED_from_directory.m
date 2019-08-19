function dataTable = read_MED_from_directory(varargin)
%% AUTHOR    : Oscar X. Guerrero-Gutierrez
%% $DATE     : 08-Feb-2019 $
%% DEVELOPED : (R2015a)
%% FILENAME  : read_MED_from_directory.m
%% Parameters
% argument parser
pArgs = inputParser;
% required parameters
% TODO: use validate functions
pArgs.addRequired('directory');
pArgs.addRequired('name_table');
pArgs.addParameter('allow_empty', false);
% aliases for key values
pArgs.addParameter('events_label', 'StartQ');
pArgs.addParameter('licks_label', 'Licks');
% parameters for analyze_MED_data.m
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
directory = args.directory;
name_table = args.name_table;
%%
dataTable = table;

argsMED = struct;
argsMED.rwd_period = args.rwd_period;
argsMED.max_IBI = args.max_IBI;
argsMED.minlag = args.minlag;
argsMED.maxlag = args.maxlag;
argsMED.resolution = args.resolution;
argsMED.L = args.L;


%%
ls_dir = dir(directory);
ls_dir = ls_dir(~ismember({ls_dir.name}, {'.', '..'}));
if isempty(ls_dir)
    if args.allow_empty
        return
    else
        error('No files found in MED directory: %s', MedPC_DATA_Directory)
    end
end
% initialize dataTable
n_files = length(ls_dir);
init_cell = cell(n_files, 1);
dataTable.filename = {ls_dir.name}';
dataTable.id = init_cell;
dataTable.Subject = init_cell;
dataTable.Group = init_cell;
dataTable.Data = init_cell;

for i = 1:length(ls_dir)
    filename = ls_dir(i).name;
    file_path = fullfile(directory, filename);
    % read MED data
    argsMED.data_MED = readMedFileData(file_path);
    argsMED.name_table = name_table;
    argsMED.events_label = args.events_label;
    argsMED.licks_label = args.licks_label;
    
    data = analyze_MED_data(argsMED);
    % save to dataTable
    % aliases
    subject = data.Parameters.Subject;
    group = data.Parameters.Group;
    date = data.Parameters.StartDate;
    % using US date format YY/DD/MM
    date_str = strcat(date(7:8), date(1:2), date(4:5));
    
    dataTable.id{i} = date_str;
    dataTable.Subject{i} = subject;
    dataTable.Group{i} = group;
    dataTable.Data{i} = data;
    
end
