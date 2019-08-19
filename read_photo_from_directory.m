function dataTable = read_photo_from_directory(varargin)
%% AUTHOR    : Oscar X. Guerrero-Gutierrez
%% $DATE     : 01-Feb-2019 $
%% DEVELOPED : (R2015a)
%% FILENAME  : analyze_photo_data.m
%% Parameters
% argument parser
pArgs = inputParser;
% parameter definition
% required parameters
% TODO: use validate functions
pArgs.addRequired('directory');
pArgs.addRequired('dataTable');
% parameters for analyze_photo_data.m
pArgs.addParameter('baseline_start', -2);
pArgs.addParameter('baseline_end', 0);
pArgs.addParameter('rwd_period', 4);
pArgs.addParameter('window_start', -3);
pArgs.addParameter('window_end', 6);
% function for peak latency
pArgs.addParameter('lat_fun', struct);
% parse arguments
pArgs.parse(varargin{:});
args = pArgs.Results;
directory = args.directory;
dataTable = args.dataTable;
lat_fun = args.lat_fun;
% args for analyze_photo_data.m
argsPhoto.baseline_start = args.baseline_start;
argsPhoto.baseline_end = args.baseline_end;
argsPhoto.rwd_period = args.rwd_period;
argsPhoto.window_start = args.window_start;
argsPhoto.window_end = args.window_end;
%% main
ls_photo_dir = dir(directory);
ls_photo_dir = ls_photo_dir(~ismember({ls_photo_dir.name}, {'.', '..'}));
if isempty(ls_photo_dir)
    error('No files found in MED path: %s', directory)
end
% TODO: check file consistency with MED data

for i = 1:length(ls_photo_dir)
    fname = ls_photo_dir(i).name;
    % subject name must be 3 character long and be at the start of filename
    % date must be at end of filename
    % "S01_<other>_YYMMDD.csv"
    subject = fname(1:3);
    date = fname(end-9:end-4);
    
    i_group = strcmp(dataTable.Subject, subject) & strcmp(dataTable.id, date);
    group = char(dataTable.Group(i_group));
    
    file_path = fullfile(directory, fname);
    T = readtable(file_path);

    % check for consistency with MED data
    assert(strcmp(dataTable.id{i}, date), 'photo data does not match MED data');

    %
    if isfield(lat_fun, group)
        argsPhoto.lat_fun = lat_fun.(group);
    elseif isfield(argsPhoto, 'lat_fun')
        argsPhoto = rmfield(argsPhoto, 'lat_fun');
    end
    
    % analyze photo data
    dataMED = dataTable.Data{i, 1};
    dataMED = analyze_photo_data(dataMED, T, argsPhoto);
    dataTable.Data{i, 1} = dataMED;
    
%     dataTable.AUC_licks{i, 1} = 0;
%     dataTable.AUC_Zscore{i, 1} = init_cell;
%     dataTable.Pearson{i} = init_cell;
%     dataTable.Z_Peak{i, 1} = init_cell;
%     dataTable.Pearson_p{i} = init_cell;
%     
    % save to dataTable
    dataTable.Data{i,1}.Values.photo = T;
    
end

end