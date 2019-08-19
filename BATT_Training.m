%% Parameters

name_table = {
    'A', 'Licks';
    'J', 'Rwd';
    'M', 'StartQ';
    'N', 'StopQ';
    'Y', 'n_Trials';
    'S', {'Water', 0.00;
          'Sac1', 0.001;
          'Sac2', 0.002;
          'Sac3', 0.003;
          'Sac4', 0.004;
          'Sac5', 0.005}
    };


Save_ppt = 0;


%% Absolute Path of MedPC DATA
MedPC_DATA_Directory = 'C:\Users\Lab36\iCloudDrive\Maestria\Lab Neurobiologia del Apetito\MedPC\Brief-access Taste Test\Batt_Training_Photometry\Other';

%Save pptx whit Rasters
Results_path = 'C:\Users\Lab36\iCloudDrive\Maestria\Lab Neurobiologia del Apetito\PHOTOMETRY_PLX\PHOTO_Results';
ppt_file = fullfile(Results_path, 'BATT_TRAINING_Raster.ppt');

%% Read and Analyze MedPC data

dataTable = read_MED_from_directory(MedPC_DATA_Directory, name_table);

%% Raster parameters

slcn_seq = {'Water' 'Sac1' 'Sac2' 'Sac3' 'Sac4' 'Sac5'};
marker_seq = {'o', '^', 's', 'd', 'v', '>'};

% color information
colorStruct.Water = 'b';  % Blue
colorStruct.Sac1 = [30,144,255]/255;  % Dodger blue
colorStruct.Sac2 = [255,215,0]/255;  % Gold
colorStruct.Sac3 = [255,165,0]/255;  % Orange
colorStruct.Sac4 = [255,105,180]/255;  % Hotpink
colorStruct.Sac5 = [220,60,60]/255;  % Crimson

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