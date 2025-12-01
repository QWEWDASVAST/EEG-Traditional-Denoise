# EEG-Traditional-Denoise
EEG Denoising with ICA, EEMD, and Wavelet Thresholding (WT)
MATLAB implementations of three traditional EEG denoising methods, supporting batch processing of PhysioNetMI, SEED, TSUBenchmark datasets.
Prerequisites
MATLAB (R2018b+ recommended)
EEGLAB Toolbox (required version: eeglab2025.0.0)
Download link: https://sccn.ucsd.edu/eeglab/download.php
Add EEGLAB to MATLAB path before running code.
Electrode Position Templates
Prepare dataset-matching template files:
biosemi_template.mat (64-channel: PhysioNetMI/TSUBenchmark)
channel_62_pos.locs (62-channel: SEED)
Place files in the script directory or update the file path in code.
