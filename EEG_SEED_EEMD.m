clear; clc;

%% Automatic batch processing settings
input_folder = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\SEED\Preprocessed_EEG';
output_folder = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\SEED\EEMD';

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

% EEMD parameter settings 
Nstd = 0.2;          % Standard deviation ratio of added noise
NE = 50;             % Number of EEMD ensembles
corr_threshold = 0.1; % IMF correlation coefficient threshold

% Define high-risk channels requiring denoising for 62-channel SEED data
% Define high-risk channels (frontal and temporal lobe regions) based on standard 62-channel system
target_channels = {'Fp1', 'Fpz', 'Fp2', 'AF3', 'AF4', 'F7', 'F8', 'T7', 'T8'}; % Core 8 channels

% Get all electrode names for 62 channels (extracted from locs)
if isstruct(locs)
    all_channel_names = {locs.labels};
else
    all_channel_names = locs;
end

fprintf('Electrode system contains %d channels\n', length(all_channel_names));

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
    
    % EEMD denoising processing
    eeg_epochs_cell = cell(1, num_datasets);
    for i = 1:num_datasets
        fprintf('  Processing dataset %d/%d...\n', i, num_datasets);
        current_eeg_data = eeg_data_cell{i};
        
        % Get indices of target channels - exactly the same as the code you provided
        [~, channels_to_denoise] = ismember(target_channels, all_channel_names);
        channels_to_denoise = channels_to_denoise(channels_to_denoise > 0);

        nTotalChans = size(current_eeg_data, 1);
        channels_to_denoise = channels_to_denoise(channels_to_denoise <= nTotalChans);

        fprintf('    Selected %d high-risk channels for EEMD denoising: ', length(channels_to_denoise));
        for j = 1:length(channels_to_denoise)
            fprintf('%s ', all_channel_names{channels_to_denoise(j)});
        end
        fprintf('\n');
        
        % Extract data of channels needing denoising
        noisy_data = current_eeg_data(channels_to_denoise, :);
        
        % Apply EEMD denoising - directly call your existing multichannel_eemd_denoise function
        denoised_data = multichannel_eemd_denoise(noisy_data, 200, Nstd, NE, corr_threshold);
        
        % Put denoised data back to original position
        current_eeg_data(channels_to_denoise, :) = denoised_data;
    
        eeg_epochs_cell{i} = current_eeg_data;
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
