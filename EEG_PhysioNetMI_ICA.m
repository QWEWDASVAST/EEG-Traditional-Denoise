clear; clc;

%% 1. Define paths and parameters
input_root = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\PhysioNetMI\files\';
output_root = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\PhysioNetMI\processing\ICA_EO_BSP1\pycache_data\';

% Load electrode position template
load('biosemi_template.mat');

% Define channels to process
EEG_channels = 1:64;

%% 2. Get all subject directories
% Use dir to get all folders starting with 'S'
subject_folders = dir(fullfile(input_root, 'S*'));
subject_folders = subject_folders([subject_folders.isdir]); 
subject_names = {subject_folders.name}; 
fprintf('Found %d subject directories\n', length(subject_names));

%% 3. Iterate through each subject directory
for sub_idx = 1:length(subject_names)
    current_subject = subject_names{sub_idx}; 
    subject_input_path = fullfile(input_root, current_subject);
    subject_output_path = fullfile(output_root, current_subject);
    
    % Create output directory if it does not exist
    if ~exist(subject_output_path, 'dir')
        mkdir(subject_output_path);
        fprintf('Created output directory: %s\n', subject_output_path);
    end
    
    % Get all EDF files in the current subject directory
    edf_files = dir(fullfile(subject_input_path, '*.edf'));
    fprintf('Processing subject %s, total %d EDF files\n', current_subject, length(edf_files));
    
    %% 4. Iterate through and process each EDF file of the current subject
    for file_idx = 1:length(edf_files)
        edf_file_name = edf_files(file_idx).name;
        edf_file_path = fullfile(subject_input_path, edf_file_name);
        
        fprintf('Start processing file: %s\n', edf_file_name);
        
        try
            % 4.1 Read EDF file
            eeg_struct = pop_biosig(edf_file_path);
            
            % 4.2 Extract data and related parameters
            fs = eeg_struct.srate; % Get original sampling rate
            
            nWant = 64;
            nChans = size(eeg_struct.data,1);
            EEG_channels = 1:min(nWant, nChans);
            EEG_data = eeg_struct.data(EEG_channels, :); % Extract data of specified channels
            
            % 4.3 Artifact removal using ICA and ICLabel
            EEG_data_cleaned = IC_label_artifact_removal(EEG_data, fs, 528, locs, EEG_channels);
            
            % 4.4 Create cleaned EEG structure
            eeg_cleaned = eeg_struct;
            eeg_cleaned.data = EEG_data_cleaned;
                     
            % 4.5 Define output path and save
            output_edf_path = fullfile(subject_output_path, edf_file_name);
            pop_writeeeg(eeg_cleaned, output_edf_path, 'TYPE', 'EDF');
            
            fprintf('Successfully processed and saved: %s\n', output_edf_path);
            
        catch ME
            % Error occurred during processing, record error message and continue processing other files
            fprintf('Error processing file %s: %s\n', edf_file_name, ME.message);
            continue; % Continue processing next file
        end
    end
end

fprintf('All files processed successfully!\n');
