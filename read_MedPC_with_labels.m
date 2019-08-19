% Made by Oscar X. Guerrero Gutierrez
% tagged variables must be in seconds
% tags must be in 0.001 multiples
function Data = read_MedPC_with_labels(filename, name_table)
    Data = readMedFileData(filename);
    for varName_label_pair = name_table'
        [varName, label] = varName_label_pair{:};
        
        % check if variable is tagged
        if iscell(label)
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
            label_cell = label;
            for label_tag_pair = label_cell'
                [label, tag] = label_tag_pair{:};
                mask = (int_tags == round(tag*1000));
                Data.Values.(label) = values(mask);
            end
        else
            raw_values = Data.Values.(varName);
            [m, n] = size(raw_values);
            values = reshape(raw_values', m*n, 1);
            % Erase empty values
            values(values==0) = [];
            values(isnan(values)) = [];
            %
            Data.Values.(label) = values;
        end
        
    end
end