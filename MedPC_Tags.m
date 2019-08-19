function [Data] = MedPC_Tags(NameFile, DIM)

% Read MPC file
Data = readMedFileData(NameFile);
% Run each DIM array
taggedVars = fieldnames(DIM);
for arrName=taggedVars
    % Extract and columnwise reshape values from ith DIM array into X
    varName = arrName{1};
    raw_values = Data.Values.(varName);
    [m, n] = size(raw_values);
    values = reshape(raw_values', m*n, 1);
    % Erase empty values
    values(values==0) = [];
    values(isnan(values)) = [];
    % Make a vector of tags, subtract it from values
    tags = mod(values, 0.01);
    int_tags = round(tags*1000);
    values = values - tags;
    
    % Generate Structure with fields named according to each event
    for label_tag_pair = DIM.(varName)'
        [label, tag] = label_tag_pair{:};
        mask = (int_tags == round(tag*1000));
        Data.Values.(label) = values(mask);
    end
    clear raw_values
    clear values
end

end %Function