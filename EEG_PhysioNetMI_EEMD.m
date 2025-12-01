clear; clc;

%% Define paths and parameters
input_root = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\PhysioNetMI\files\';
output_root = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\PhysioNetMI\processing\Phy_EEMD\pycache_data\';

% Load electrode position template
load('biosemi_template.mat');

PHYSIONETMI_CHANNEL_LIST = {'Fc5.', 'Fc3.', 'Fc1.', 'Fcz.', 'Fc2.', 'Fc4.', 'Fc6.', 'C5..', 'C3..', 'C1..', 'Cz..', 'C2..', 'C4..', 'C6..', 'Cp5.', 'Cp3.', 'Cp1.', 'Cpz.', 'Cp2.', 'Cp4.', 'Cp6.', 'Fp1.', 'Fpz.', 'Fp2.', 'Af7.', 'Af3.', 'Afz.', 'Af4.', 'Af8.', 'F7..', 'F5..', 'F3..', 'F1..', 'Fz..', 'F2..', 'F4..', 'F6..', 'F8..', 'Ft7.', 'Ft8.', 'T7..', 'T8..', 'T9..', 'T10.', 'Tp7.', 'Tp8.', 'P7..', 'P5..', 'P3..', 'P1..', 'Pz..', 'P2..', 'P4..', 'P6..', 'P8..', 'Po7.', 'Po3.', 'Poz.', 'Po4.', 'Po8.', 'O1..', 'Oz..', 'O2..', 'Iz..'};
% EEMD parameter settings
Nstd = 0.2;          % Standard deviation ratio of added noise
NE = 50;             % Number of EEMD ensembles
corr_threshold = 0.1; % IMF correlation coefficient threshold

% Define high-risk channels requiring denoising (frontal and temporal lobe regions)
target_channels = {'Fp1.', 'Fpz.', 'Fp2.', 'Af7.', 'Af8.', 'F7..', 'F8..', 'Iz..'}; % Core 8 channels
% target_channels = {'Fp1', 'Fpz', 'Fp2', 'Af7', 'Af8', 'F7', 'F8', 'Ft7', 'Ft8', 'Iz'}; % Extended 10 channels

% Define channels to process
EEG_channels = 1:64;

%% Get all subject directories
% Use dir to get all folders starting with 'S'
subject_folders = dir(fullfile(input_root, 'S*'));
subject_folders = subject_folders([subject_folders.isdir]); % Keep only directories, exclude files
subject_names = {subject_folders.name}; % Get cell array of directory names
% subject_names = {'S106'}; % Get cell array of directory names (specific subject for testing)
fprintf('Found %d subject directories\n', length(subject_names));

%%  Iterate through each subject directory
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
    
    %% Iterate through and process each EDF file of the current subject
    for file_idx = 1:length(edf_files)
        edf_file_name = edf_files(file_idx).name;
        edf_file_path = fullfile(subject_input_path, edf_file_name);
        
        fprintf('Start processing file: %s\n', edf_file_name);
        
        try
            eeg_struct = pop_biosig(edf_file_path);
            
            %Extract data and related parameters
            fs = eeg_struct.srate; % Get original sampling rate
            
            nWant = 64;
            nChans = size(eeg_struct.data, 1);
            EEG_channels = 1:min(nWant, nChans);
            EEG_data = eeg_struct.data(EEG_channels, :); % Extract data of specified channels
            
            % Artifact removal using EEMD
            [~, channels_to_denoise] = ismember(target_channels, PHYSIONETMI_CHANNEL_LIST);
            channels_to_denoise = channels_to_denoise(channels_to_denoise > 0);
            
            nTotalChans = size(EEG_data, 1);
            channels_to_denoise = channels_to_denoise(channels_to_denoise <= nTotalChans);
            
            fprintf('Selected %d high-risk channels for EEMD denoising:\n', length(channels_to_denoise));
            for i = 1:length(channels_to_denoise)
                fprintf('%s ', PHYSIONETMI_CHANNEL_LIST{channels_to_denoise(i)});
            end
            fprintf('\n');
            
            %Perform EEMD denoising only on high-risk channels
            EEG_data_cleaned = EEG_data; % Initialize output data
            
            if ~isempty(channels_to_denoise)
                fprintf('Starting EEMD denoising...\n');
                % Extract data of channels needing denoising
                noisy_data = EEG_data(channels_to_denoise, :);
                % Apply EEMD denoising
                denoised_data = multichannel_eemd_denoise(noisy_data, fs, Nstd, NE, corr_threshold);
                % Put denoised data back to original position
                EEG_data_cleaned(channels_to_denoise, :) = denoised_data;
            end
            
            clean_channels = setdiff(1:size(EEG_data, 1), channels_to_denoise);
            fprintf('The remaining %d channels retain original signals.\n', length(clean_channels));
            
%             EEG_data_cleaned = multichannel_eemd_denoise(EEG_data, fs, target_fs, Nstd, NE, corr_threshold);
            
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
