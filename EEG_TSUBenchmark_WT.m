clear; clc;

%% Prepare file paths
data_path = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\TSUBenchmark\Raw';
new_save_path = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\TSUBenchmark\processing\WT';

% Ensure output directory exists
if ~isfolder(new_save_path)
    mkdir(new_save_path);
end

%% Get list of all .mat files
% Use dir function to get all .mat files in the directory
mat_files = dir(fullfile(data_path, '*.mat'));
fprintf('Found %d MAT files to process\n', length(mat_files));

fs = 250;

% Wavelet threshold denoising parameters
wname = 'db4';        % Wavelet name, db4 wavelet is commonly used in EEG processing
level = 3;            % Wavelet decomposition level

%% Load electrode position template
load('biosemi_template.mat');          % Variables locs / EEG_channels already exist
topoplot([], locs, 'electrodes','ptslabels','plotdisk','on');

EEG_channels = 1:64;

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
        
        %% Process each epoch individually with wavelet threshold denoising
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
                
                % Perform wavelet threshold denoising on current epoch
                EEG_clean_epoch = WT_artifact_removal(current_epoch, EEG_channels, wname, level);
                
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

%% Wavelet threshold denoising function
function EEG_data_cleaned = WT_artifact_removal(EEG_data, EEG_channels, wname, level)
    % Perform artifact removal using wavelet thresholding (standard wavelet threshold denoising method)
    fprintf('    Performing wavelet threshold denoising...\n');
    
    % Extract data of specified channels
    EEG_data_selected = EEG_data(EEG_channels, :);
    EEG_data_cleaned = zeros(size(EEG_data_selected));
    
    % Perform wavelet threshold denoising for each channel individually
    for ch = 1:size(EEG_data_selected, 1)
        channel_data = EEG_data_selected(ch, :);
        N = length(channel_data);
        
        % Use standard wavelet threshold denoising method
        % Wavelet decomposition
        [c, l] = wavedec(channel_data, level, wname);
        
        % Apply soft thresholding to each detail coefficient
        for i = 1:level
            start_index = sum(l(1:i)) + 1;
            end_index = sum(l(1:i+1));
            
            % Use universal threshold formula
            thr = 0.5 * sqrt(2 * log(N)); % Threshold: 0.5*sqrt(2*log(N))
            
            % Apply soft thresholding
            c(start_index:end_index) = wthresh(c(start_index:end_index), 's', thr); 
        end
        
        % Reconstruct denoised signal using wavelet
        cleaned_channel = waverec(c, l, wname);
        
        EEG_data_cleaned(ch, :) = cleaned_channel;
    end

end
