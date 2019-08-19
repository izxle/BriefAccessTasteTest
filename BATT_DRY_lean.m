%% Parameters

name_table = {
    'A', 'Licks';
    'J', 'Rwd';
    'M', 'StartQ';
    'N', 'StopQ';
    'Y', 'n_Trials';
    'S', {'Water', 0.00;
          'Dry', 0.001;
          'Sac1', 0.002;
          'Sac2', 0.003}
    };

Save_ppt = 0;


%% Absolute Path of MedPC DATA
MedPC_DATA_Directory = 'C:\Users\Lab36\iCloudDrive\Maestria\Lab Neurobiologia del Apetito\\MedPC\Brief-access Taste Test\Batt_DryLicks_Photometry';

PHOTO_DATA_Directory = 'C:\Users\Lab36\iCloudDrive\Maestria\Lab Neurobiologia del Apetito\PHOTOMETRY_PLX\Mona PHOTO\BATT_Dry_Photo';

%Save pptx whit Rasters
Results_path = 'C:\Users\Mona\iCloudDrive\Maestria\Lab Neurobiologia del Apetito\PHOTOMETRY_PLX\PHOTO_Results';
ppt_file = fullfile(Results_path, 'BATT_DRY_PHOTOMETRY_Raster_2.ppt');

%% Read and Analyze MedPC data

dataTable = read_MED_from_directory(MedPC_DATA_Directory, name_table);

%% Read Photometry data
photoArgs.baseline_start = -3;
photoArgs.baseline_end = -1;
photoArgs.window_start = -3;
photoArgs.window_end = 6;

photo_x_time = linspace(-3, 6, 270);

photoArgs.lat_fun = struct('Vgat_GCaMP6', 'max', 'Vgat_Vector', 'min', 'Vglut2_GCaMP6', 'min');
dataTable = read_photo_from_directory(PHOTO_DATA_Directory, dataTable, photoArgs);

%% Quantilify data
quantArgs = struct;
quantArgs.p = [0.25 0.75];

dataTable = quantilify_data(dataTable, quantArgs);

%% Quiantilify Photometry data
args = struct;
args.photo_x_time = photo_x_time;
args.photo_Zscore_y_lim = [-2 4];

colors.q1 = 'b';  % Blue
colors.q2 = [220,60,60]/255;  % Crimson

photo_by_quantile(dataTable, args);


%% Behaviour data Stratifies ANOVA & Pos hoc Fisher
% slcn_seq = {'Dry' 'Water' 'Sac1' 'Sac2'};
% 
% palatables = BATT_stratifydata(dataTable, slcn_seq);
% palTable = dataTable(palatables, :);
% 
% notpalTable = dataTable(~palatables, :);

%% Parameters for BATT Descriptors
slcn_seq = {'Dry' 'Water' 'Sac1' 'Sac2'};
marker_seq = {'o', '^', 's', 'd', 'v', '>', 'h', '.', '+'};

% color information
colorStruct.Water = 'b';  % Blue
colorStruct.Dry = [105,105,105]/255;  % Dimgray
colorStruct.Sac1 = [255,215,0]/255;  % Gold
colorStruct.Sac2 = [255,105,180]/255;  % Hotpink

% Group color information, name the group
colorGroup.Vgat_GCaMP6 = [57, 83, 164]/255;  % Medium Purple
colorGroup.Vgat_Vector = [105, 105, 105]/255; % Dimgray
colorGroup.Vglut2_GCaMP6 = [222, 45, 38]/255;  % Burly Wood
colorGroup.Vgat_eGFP = [128, 128, 128]/255;

descriptorsArgs.color_slcn = colorStruct;
descriptorsArgs.color_group = colorGroup;
descriptorsArgs.marker_size = 110;
descriptorsArgs.marker_seq = marker_seq;
descriptorsArgs.photo_x_time = linspace(-3, 6, 270);

descriptorsArgs.bar_x_label = {'Empty', 'Water', '3','18'};

descriptorsArgs.photo_Zscore_y_lim = [-2 4];

%%
BATT_Descriptors(dataTable, slcn_seq, descriptorsArgs);

%%
licks_outArgs.color = [211,211,211]/255;
licks_outArgs.line_style_seq = {'-', '--'};

photo_licks_after(dataTable, licks_outArgs);

%%  plot raster for all files
n_files = size(dataTable, 1);
for i_file = 1:n_files
    data = dataTable.Data{i_file};
    % build arguments for raster
    date = data.Parameters.StartDate;
    subject = data.Parameters.Subject;
    fig_title = strcat(subject, '-', date);
    rasterArgs.fig_title = fig_title; 
    rasterArgs.rwd_period = 4;
    rasterArgs.slcn_seq = slcn_seq;  
    rasterArgs.colorStruct = colorStruct;
    rasterArgs.color_rwd = [180,67,212]/255;  % "Lila"
    rasterArgs.L = 5;
    %
%     Raster_BATT_photo(data, rasterArgs);

    %Save Raster/PSTH on .ppt   
    if Save_ppt
       saveppt2(ppt_file);
       close all
    end
end