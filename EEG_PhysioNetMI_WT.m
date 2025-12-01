clear; clc;
addpath(genpath('E:\workspace\Pycharmworkspace\EEGdenoiseNet-master\code\preprocessing\EEG\waveletTh_wp'))

%% Define paths and parameters
input_root = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\PhysioNetMI\files\';
output_root = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\PhysioNetMI\processing\WT\pycache_data\';

% Load electrode position template
load('biosemi_template.mat'); 

% Define channels to process
EEG_channels = 1:64;
% Wavelet threshold denoising parameters
wname = 'db4';        % Wavelet name, db4 wavelet is commonly used in EEG processing
level = 3;            % Wavelet decomposition level


%% Get all subject directories
% Use dir to get all folders starting with 'S'
subject_folders = dir(fullfile(input_root, 'S*'));
subject_folders = subject_folders([subject_folders.isdir]); 
subject_names = {subject_folders.name}; % Get cell array of directory names
% subject_names = {'S106'}; % Get cell array of directory names (specific subject for testing)
fprintf('Found %d subject directories\n', length(subject_names));

%% 3. Iterate through each subject directory
for sub_idx = 1:length(subject_names)
    current_subject = subject_names{sub_idx}; 
    subject_input_path = fullfile(input_root, current_subject);
    subject_output_path = fullfile(output_root, current_subject);
    
    % Create output directory
    if ~exist(subject_output_path, 'dir')
        mkdir(subject_output_path);
        fprintf('Created output directory: %s\n', subject_output_path);
    end
    
    % Get all EDF files in the current subject directory
    edf_files = dir(fullfile(subject_input_path, '*.edf'));
    fprintf('Processing subject %s, total %d EDF files\n', current_subject, length(edf_files));
    
    %% Iterate through and process each EDF file of the current subject
    for file_idx = 1:length(edf_files)
        edf_file_name = edf_files(file_idx).name;
        edf_file_path = fullfile(subject_input_path, edf_file_name);
        
        fprintf('Start processing file: %s\n', edf_file_name);
        
        try
            % Read EDF file
            eeg_struct = pop_biosig(edf_file_path);
            
            % Extract data and related parameters
            fs = eeg_struct.srate; % Get original sampling rate
            
            nWant = 64;
            nChans = size(eeg_struct.data, 1);
            EEG_channels = 1:min(nWant, nChans);
            EEG_data = eeg_struct.data(EEG_channels, :); % Extract data of specified channels
            
            % Artifact removal using wavelet thresholding (standard wavelet threshold denoising method)
            fprintf('Performing wavelet threshold denoising...\n');
            EEG_data_cleaned = zeros(size(EEG_data));
            
            % Perform wavelet threshold denoising for each channel individually
            for ch = 1:size(EEG_data, 1)
                channel_data = EEG_data(ch, :);
                N = length(channel_data);
                
                % Use standard wavelet threshold denoising method
                % Wavelet decomposition
                [c, l] = wavedec(channel_data, level, wname);
                
                % Apply soft thresholding to each detail coefficient
                for i = 1:level
                    start_index = l(i) + 1;
                    end_index = l(i + 1);
                    
                    % Use universal threshold formula
                    thr = 0.5 * sqrt(2 * log(N)); % Threshold: 0.5*sqrt(2*log(N))
                    
                    % Apply soft thresholding
                    c(start_index:end_index) = wthresh(c(start_index:end_index), 's', thr); 
                end
                
                % Reconstruct denoised signal using wavelet
                cleaned_channel = waverec(c, l, wname);
                
                EEG_data_cleaned(ch, :) = cleaned_channel;
            end
            
            % Create cleaned EEG structure
            eeg_cleaned = eeg_struct; % Copy original structure
            eeg_cleaned.data = EEG_data_cleaned; % Replace with cleaned data
           
            
            % Define output path and save
            output_edf_path = fullfile(subject_output_path, edf_file_name);
            pop_writeeeg(eeg_cleaned, output_edf_path, 'TYPE', 'EDF');
            
            fprintf('Successfully processed and saved: %s\n', output_edf_path);
            
        catch ME
            % Error occurred during processing, record error message and continue with other files
            fprintf('Error processing file %s: %s\n', edf_file_name, ME.message);
            continue; % Continue processing next file
        end
    end
end

fprintf('All files processed successfully!\n');
