clear; clc;

%% Prepare file paths
data_path = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\TSUBenchmark\Raw';
new_save_path = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\TSUBenchmark\processing\EEMD';

% Ensure output directory exists
if ~isfolder(new_save_path)
    mkdir(new_save_path);
end

%% Get list of all .mat files
% Use dir function to get all .mat files in the directory
mat_files = dir(fullfile(data_path, '*.mat'));
fprintf('Found %d MAT files to process\n', length(mat_files));

fs = 250;

%% Load electrode position template
load('biosemi_template.mat');          % Variables locs / EEG_channels already exist
topoplot([], locs, 'electrodes','ptslabels','plotdisk','on');
PHYSIONETMI_CHANNEL_LIST = {'Fc5.', 'Fc3.', 'Fc1.', 'Fcz.', 'Fc2.', 'Fc4.', 'Fc6.', 'C5..', 'C3..', 'C1..', 'Cz..', 'C2..', 'C4..', 'C6..', 'Cp5.', 'Cp3.', 'Cp1.', 'Cpz.', 'Cp2.', 'Cp4.', 'Cp6.', 'Fp1.', 'Fpz.', 'Fp2.', 'Af7.', 'Af3.', 'Afz.', 'Af4.', 'Af8.', 'F7..', 'F5..', 'F3..', 'F1..', 'Fz..', 'F2..', 'F4..', 'F6..', 'F8..', 'Ft7.', 'Ft8.', 'T7..', 'T8..', 'T9..', 'T10.', 'Tp7.', 'Tp8.', 'P7..', 'P5..', 'P3..', 'P1..', 'Pz..', 'P2..', 'P4..', 'P6..', 'P8..', 'Po7.', 'Po3.', 'Poz.', 'Po4.', 'Po8.', 'O1..', 'Oz..', 'O2..', 'Iz..'};
EEG_channels = 1:64;

% Define high-risk channels requiring denoising (frontal and temporal lobe regions)
target_channels = {'Fp1.', 'Fpz.', 'Fp2.', 'Af7.', 'Af8.', 'F7..', 'F8..', 'Iz..'};

% EEMD parameter settings
Nstd = 0.2;          % Standard deviation ratio of added noise
NE = 50;             % Number of EEMD ensembles
corr_threshold = 0.1; % IMF correlation coefficient threshold

%% Process each file in loop
for file_idx = 1:length(mat_files)
    try
        % Get current file name
        data_file = mat_files(file_idx).name;
        fid = fullfile(data_path, data_file);
        
        fprintf('\nProcessing file %d/%d: %s\n', file_idx, length(mat_files), data_file);
        
        % Load data
        load(fid);
        
        fs = 250;
        
        %% Parameter initialization
        [nChan, nTime, nFreq, nBlock] = size(data);
        nEpoch = nFreq * nBlock;               % Total number of epochs
        
        % Initialize 4D array for denoised data
        data_denoised = zeros(nChan, nTime, nFreq, nBlock);
        
        %% Process each epoch individually with EEMD artifact removal
        fprintf('  Starting individual epoch processing, total %d epochs...\n', nEpoch);
        
        % Loop through each epoch
        epoch_count = 0;
        for freq_idx = 1:nFreq
            for block_idx = 1:nBlock
                epoch_count = epoch_count + 1;
                fprintf('  Processing epoch %d/%d (frequency %d/%d, block %d/%d)...\n', ...
                    epoch_count, nEpoch, freq_idx, nFreq, block_idx, nBlock);
                
                % Extract data of current epoch (64, 1500)
                current_epoch = data(:, :, freq_idx, block_idx);
                
                % 4.3 Artifact removal using EEMD
                [~, channels_to_denoise] = ismember(target_channels, PHYSIONETMI_CHANNEL_LIST);
                channels_to_denoise = channels_to_denoise(channels_to_denoise > 0);

                nTotalChans = size(current_epoch, 1);
                channels_to_denoise = channels_to_denoise(channels_to_denoise <= nTotalChans);

                fprintf('Selected %d high-risk channels for EEMD denoising:\n', length(channels_to_denoise));
                
                EEG_clean_epoch = current_epoch;
                
                % Extract data of channels needing denoising
                noisy_data = EEG_clean_epoch(channels_to_denoise, :);
                % Apply EEMD denoising
                denoised_data = multichannel_eemd_denoise(noisy_data, fs, Nstd, NE, corr_threshold);
                % Put denoised data back to original position
                EEG_clean_epoch(channels_to_denoise, :) = denoised_data;
               
                % Verify dimensions of processed data
                [clean_nChan, clean_nTime] = size(EEG_clean_epoch);
                if clean_nChan ~= nChan || clean_nTime ~= nTime
                    fprintf('    Warning: Epoch dimension changed after processing (%d×%d -> %d×%d), performing automatic adjustment\n', ...
                        nChan, nTime, clean_nChan, clean_nTime);
                    % Adjust dimensions or use original data if mismatch occurs
                    if clean_nChan == nChan && clean_nTime <= nTime
                        % Zero-padding if channel count matches but time points are reduced
                        temp_data = zeros(nChan, nTime);
                        temp_data(:, 1:min(nTime, clean_nTime)) = EEG_clean_epoch(:, 1:min(nTime, clean_nTime));
                        EEG_clean_epoch = temp_data;
                    else
                        fprintf('    Error: Unable to adjust dimensions, using original epoch data\n');
                        EEG_clean_epoch = current_epoch;
                    end
                end
                
                % Store processed epoch data back to 4D array
                data_denoised(:, :, freq_idx, block_idx) = EEG_clean_epoch;
            end
        end
        
        fprintf('  All epochs processed successfully!\n');
        
        %% Save processed data
        full_file_path = fullfile(new_save_path, data_file);
        % Load all variables from original file
        file_vars = load(fid);
        vars_list = fieldnames(file_vars);
        
        % Update the variable to be modified
        file_vars.data = data_denoised;  % Replace with denoised data
        
        % Save all variables (including other potential variables in original file)
        save(full_file_path, '-struct', 'file_vars');
        fprintf('   All variable structures from original file are preserved, and data has been updated\n');
        
    catch ME
        fprintf('Error: An error occurred while processing file %s: %s\n', data_file, ME.message);
        fprintf('Error stack trace:\n');
        for k = 1:length(ME.stack)
            fprintf('  %s (%d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        continue; % Skip current file and proceed to next one
    end
end

fprintf('\nBatch processing completed! Processed %d/%d files\n', file_idx, length(mat_files));
