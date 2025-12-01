% This function references the preprocessing code from the EEGdenoiseNet dataset
function [EEG_data_processed] = IC_label_artifact_removal(EEG_data, fs, order, locs, EEG_channels)
% Add paths for required libraries and FastICA implementation, both sourced from the EEGdenoiseNet project's preprocessing module
% Project GitHub URL: https://github.com/ncclabsustech/EEGdenoiseNet/tree/master/code/preprocessing

% Add path for basic auxiliary libraries required for EEG preprocessing
% (Corresponds to the content in the 'preprocessing/lib' directory of the GitHub repository)
% This directory contains various helper functions necessary for EEG signal preprocessing workflows
addpath(genpath('E:\workspace\Pycharmworkspace\EEGdenoiseNet-master\code\preprocessing\lib'));

% Add path for FastICA algorithm implementation
% (Corresponds to the content in the 'preprocessing/EEG/fastica' directory of the GitHub repository)
% This directory includes the MATLAB implementation of FastICA (Fast Independent Component Analysis),
% which is used for independent component decomposition of EEG signals
addpath(genpath('E:\workspace\Pycharmworkspace\EEGdenoiseNet-master\code\preprocessing\EEG\fastica'));

%% 2. Re-reference and bandpass filtering
low_cutoff = 1;    
% Automatically set high-pass cutoff frequency based on sampling rate: 79 Hz for 160 Hz, 63 Hz for 128 Hz, 80 Hz otherwise
if fs == 160
    high_cutoff = 79;      
elseif fs == 128
    high_cutoff = 63;     
else
    high_cutoff = 80;     
end
% Bandpass filter (1 - high_cutoff Hz) to remove very low-frequency drift and high-frequency noise (e.g., EMG)
EEG_data = filter_data(EEG_data, fs, low_cutoff, high_cutoff, order);

%% 4. ICA decomposition using FastICA
% 3.1 PCA whitening
[dewhitening, score, latent] = pca(EEG_data', 'Centered', 'off');
whitening = inv(dewhitening);
PC = whitening * EEG_data;

% Calculate variance contribution rate of principal components, 
% determine the number of retained PCs (n_pcs) based on threshold (0.01% of average variance)
% This is used for dimensionality reduction to reduce computation and avoid overfitting
latent = 100 * latent / sum(latent);
n_pcs = sum(latent > 0.01 * mean(latent));
retained_variance = sum(latent(1:n_pcs));

% Perform FastICA decomposition
[IC, mixing, unmixing] = fastica(PC(1:n_pcs, :), 'approach', 'defl', 'g', 'tanh', 'maxNumIterations', 500); 
IC = unmixing * PC(1:n_pcs, :);
mixing_matrix = dewhitening(:, 1:n_pcs) * mixing;
unmixing_matrix = unmixing * whitening(1:n_pcs, :);

%% 4. Use ICLabel for artifact classification
% Prepare EEGLAB standard data structure for ICLabel
data_len = size(EEG_data, 2);
time_len = data_len / fs;
time_axis = linspace(0, time_len, data_len);
EEG_struct.times = time_axis;
EEG_struct.data = EEG_data;
EEG_struct.chanlocs = locs;
EEG_struct.srate = fs;
EEG_struct.trials = 1;
EEG_struct.pnts = data_len;
EEG_struct.icawinv = mixing_matrix;
EEG_struct.icaweights = unmixing_matrix;
EEG_struct.icaact = unmixing_matrix * EEG_data;
EEG_struct.icachansind = EEG_channels;
EEG_struct.ref = 'averef';

% Perform ICLabel classification
EEG_struct = iclabel(EEG_struct);
class_prob = EEG_struct.etc.ic_classification.ICLabel.classifications;
class_type = EEG_struct.etc.ic_classification.ICLabel.classes;

%% 7. Reconstruct cleaned EEG data
ic_num = size(class_prob, 1);
for iter_ic = 1:ic_num
    [max_prob_val(1, iter_ic), class_index(1, iter_ic)] = max(class_prob(iter_ic, :));
end

% Categorize ICs based on classification results
class_label.brain_ic_index = sort(find(class_index == 1));
class_label.muscle_ic_index = sort(find(class_index == 2));
class_label.eye_ic_index = sort(find(class_index == 3));
class_label.heart_ic_index = sort(find(class_index == 4));
class_label.line_noise_ic_index = sort(find(class_index == 5));
class_label.channel_noise_ic_index = sort(find(class_index == 6));
class_label.other_noise_ic_index = sort(find(class_index == 7)); % "Other" category - consider whether to remove based on specific needs

brain_ic_num = length(class_label.brain_ic_index);
if (brain_ic_num == 0)
    warning('No brain ICs recovered! Please check the data or parameters.');
end

% Identify bad ICs (artifacts) to remove
bad_ic_index = unique(sort([class_label.muscle_ic_index, class_label.eye_ic_index, ...
    class_label.heart_ic_index, class_label.line_noise_ic_index, ...
    class_label.channel_noise_ic_index])); % Uncomment class_label.other_noise_ic_index if needed

% Reconstruct cleaned EEG by removing bad ICs
EEG_data_processed = EEG_data - mixing_matrix(:, bad_ic_index) * IC(bad_ic_index, :);
end
