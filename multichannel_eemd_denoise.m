function [denoised_data] = multichannel_eemd_denoise(EEG_data, fs, Nstd, NE, corr_threshold)
    % Function: Perform EEMD denoising on multichannel EEG data
    % Inputs:
    %   EEG_data: Multichannel EEG data (channels x samples)
    %   fs: Original sampling rate
    %   Nstd: Standard deviation ratio of added noise
    %   NE: Number of EEMD ensembles
    %   corr_threshold: IMF correlation coefficient threshold
    % Output:
    %   denoised_data: Denoised EEG data
    
    [num_channels, num_samples] = size(EEG_data);
    denoised_data = zeros(num_channels, num_samples);
    
    fprintf('Starting multichannel EEMD denoising with correlation coefficient threshold: %.3f\n', corr_threshold);
    fprintf('Parameter settings - Nstd: %.2f, NE: %d\n', Nstd, NE);
    
    % Check parallel pool status; start if not running
    if isempty(gcp('nocreate'))
        try
            parpool('local'); % Use default number of workers
            fprintf('Parallel pool started successfully\n');
        catch ME
            fprintf('Failed to start parallel pool, will use serial processing: %s\n', ME.message);
        end
    else
        pool = gcp;
        fprintf('Parallel pool already exists with %d workers available\n', pool.NumWorkers);
    end
    
    % Process each channel in parallel using parfor
    parfor ch = 1:num_channels
        channel_start_time = tic;
        
        try
            % Extract current channel data
            channel_data = EEG_data(ch, :);
            
            % Perform EEMD decomposition
            allmode = eemd(channel_data, Nstd, NE);
            
            % Get number of components
            num_components = size(allmode, 2);
            num_imfs = num_components - 2; % Subtract original signal and residual
            
            if num_imfs <= 0
                warning('Channel %d: Failed to decompose valid IMFs, using original signal.', ch);
                denoised_channel = channel_data;
            else
                % Calculate correlation coefficient between each IMF and original signal
                correlations = zeros(1, num_imfs);
                original_signal = allmode(:, 1);
                
                for i = 1:num_imfs
                    imf_signal = allmode(:, i+1);
                    corr_matrix = corrcoef(original_signal, imf_signal);
                    correlations(i) = abs(corr_matrix(1, 2));
                end
                
                % Select IMFs to retain based on correlation coefficient threshold
                imfs_to_keep = find(correlations > corr_threshold);
                
                if isempty(imfs_to_keep)
                    % If all IMF correlations are too low, retain top few most correlated ones
                    [~, sorted_indices] = sort(correlations, 'descend');
                    imfs_to_keep = sorted_indices(1:min(3, num_imfs));
                end
                
                % Reconstruct denoised signal
                imf_indices = imfs_to_keep + 1;
                residual_index = size(allmode, 2);
                denoised_channel = sum(allmode(:, [imf_indices, residual_index]), 2)';
                
                % Log processing information (fprintf is safe in parfor)
                fprintf('Channel %d: Retained %d/%d IMFs, maximum correlation coefficient: %.3f\n', ...
                    ch, length(imfs_to_keep), num_imfs, max(correlations));
            end
            
            % Store denoised channel data
            denoised_data(ch, :) = denoised_channel;
            
            channel_time = toc(channel_start_time);
            fprintf('Channel %d processing completed, processing time: %.2f seconds\n', ch, channel_time);
            
        catch ME
            % If processing fails, use original signal and log error
            fprintf('Channel %d processing failed: %s, using original signal\n', ch, ME.message);
            denoised_data(ch, :) = EEG_data(ch, :);
        end
    end
    
    fprintf('Multichannel EEMD denoising completed\n');
end
