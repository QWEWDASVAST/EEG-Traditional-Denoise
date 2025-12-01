clear; clc;

%% Automatic batch processing settings
input_folder = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\SEED\Preprocessed_EEG';
output_folder = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\SEED\ICA_EO_BSP';

% Create output folder
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    fprintf('Created output directory: %s\n', output_folder);
end

% Get list of all .mat files to process (exclude label.mat)
all_files = dir(fullfile(input_folder, '*.mat'));
file_list = {};
for i = 1:length(all_files)
    if ~strcmp(all_files(i).name, 'label.mat')
        file_list{end+1} = all_files(i).name;
    end
end

fprintf('Found %d files to process\n', length(file_list));

% Set EEG channels
EEG_channels = 1:62;

% Prepare electrode positions
locs = readlocs('channel_62_pos.locs');
if size(locs, 1) > 62
    locs = locs(1:62, :);
    fprintf('Electrode positions adjusted to 62 channels\n');
end

%% Batch process all files
for file_idx = 1:length(file_list)
    current_file = file_list{file_idx};
    fprintf('\nProcessing file %d/%d: %s\n', file_idx, length(file_list), current_file);
    
    % Load current data file
    input_path = fullfile(input_folder, current_file);
    loaded_data = load(input_path);
    
    % Automatically detect variable name pattern
    field_names = fieldnames(loaded_data);
    fprintf('  Fields in file: %s\n', strjoin(field_names, ', '));
    
    % Detect variable name pattern using regular expression
    eeg_pattern = '^([a-zA-Z]+)_eeg\d+$';
    eeg_fields = {};
    
    for i = 1:length(field_names)
        if ~isempty(regexp(field_names{i}, eeg_pattern, 'once'))
            eeg_fields{end+1} = field_names{i};
        end
    end
    
    if isempty(eeg_fields)
        error('No variables matching *_eeg* pattern found in file %s', current_file);
    end
    
    % Extract prefix and sort fields
    prefix = regexp(eeg_fields{1}, '^([a-zA-Z]+)_eeg', 'tokens');
    prefix = prefix{1}{1};
    
    % Sort EEG fields by number
    eeg_numbers = cellfun(@(x) str2double(regexp(x, '\d+$', 'match')), eeg_fields);
    [~, sort_idx] = sort(eeg_numbers);
    eeg_fields_sorted = eeg_fields(sort_idx);
    
    num_datasets = length(eeg_fields_sorted);
    fprintf('  Detected variable pattern: %s_eeg*, total %d datasets\n', prefix, num_datasets);
    
    % Extract data into cell array
    eeg_data_cell = cell(1, num_datasets);
    for i = 1:num_datasets
        var_name = eeg_fields_sorted{i};
        eeg_data_cell{i} = loaded_data.(var_name);
        fprintf('  Loaded %s: size %dx%d\n', var_name, size(eeg_data_cell{i}, 1), size(eeg_data_cell{i}, 2));
    end
    
    % Process all datasets using ICA and ICLabel method
    eeg_epochs_cell = cell(1, num_datasets);
    for i = 1:num_datasets
        fprintf('  Processing dataset %d/%d...\n', i, num_datasets);
        current_eeg_data = eeg_data_cell{i};
        eeg_epochs_cell{i} = IC_label_artifact_removal(current_eeg_data, 200, 660, locs, EEG_channels);
    end
    
    % Save processed data, maintaining the same structure and filename as the original
    output_file = fullfile(output_folder, current_file);

    % Create structure to store all processed data
    output_data = struct();
    for i = 1:num_datasets
        var_name = eeg_fields_sorted{i};
        output_data.(var_name) = eeg_epochs_cell{i};
    end

    % Save using save function, avoiding v7.3 format
    save(output_file, '-struct', 'output_data');

    
    fprintf('File successfully saved: %s\n', output_file);
    fprintf('Saved variables: %s\n', strjoin(eeg_fields_sorted, ', '));
end

fprintf('\nBatch processing completed! Total %d files processed\n', length(file_list));
