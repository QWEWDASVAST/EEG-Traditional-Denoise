clear; clc;
target_fs = 500; % Hz 目标采样率
% raw_fs = 512; % 原始采样率将从每个文件中读取

%% 1. 定义路径和参数
input_root = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\PhysioNetMI\files\';
output_root = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\PhysioNetMI\processing\Phy_EEMD\pycache_data\';

% 加载电极位置模板（假设所有文件共用同一套电极位置）
load('biosemi_template.mat'); % 确保此文件在MATLAB路径中

PHYSIONETMI_CHANNEL_LIST = {'Fc5.', 'Fc3.', 'Fc1.', 'Fcz.', 'Fc2.', 'Fc4.', 'Fc6.', 'C5..', 'C3..', 'C1..', 'Cz..', 'C2..', 'C4..', 'C6..', 'Cp5.', 'Cp3.', 'Cp1.', 'Cpz.', 'Cp2.', 'Cp4.', 'Cp6.', 'Fp1.', 'Fpz.', 'Fp2.', 'Af7.', 'Af3.', 'Afz.', 'Af4.', 'Af8.', 'F7..', 'F5..', 'F3..', 'F1..', 'Fz..', 'F2..', 'F4..', 'F6..', 'F8..', 'Ft7.', 'Ft8.', 'T7..', 'T8..', 'T9..', 'T10.', 'Tp7.', 'Tp8.', 'P7..', 'P5..', 'P3..', 'P1..', 'Pz..', 'P2..', 'P4..', 'P6..', 'P8..', 'Po7.', 'Po3.', 'Poz.', 'Po4.', 'Po8.', 'O1..', 'Oz..', 'O2..', 'Iz..'};
% EEMD参数设置
Nstd = 0.2;          % 添加噪声的标准差比率
NE = 50;             % EEMD集成次数
corr_threshold = 0.1; % IMF相关系数阈值

% 定义需要去噪的高风险通道（前额和颞叶区域）
target_channels = {'Fp1.', 'Fpz.', 'Fp2.', 'Af7.', 'Af8.', 'F7..', 'F8..', 'Iz..'}; % 核心8通道
% target_channels = {'Fp1', 'Fpz', 'Fp2', 'Af7', 'Af8', 'F7', 'F8', 'Ft7', 'Ft8', 'Iz'}; % 扩展10通道

% 定义要处理的通道
EEG_channels = 1:64;

%% 2. 获取所有被试目录
% 使用dir获取所有以'S'开头的文件夹（假设被试目录名为100,106）
subject_folders = dir(fullfile(input_root, 'S*'));
subject_folders = subject_folders([subject_folders.isdir]); % 只保留目录，排除文件
subject_names = {subject_folders.name}; % 获取目录名称单元格数组
% subject_names = {'S106'}; % 获取目录名称单元格数组
fprintf('找到 %d 个被试目录\n', length(subject_names));

%% 3. 遍历每个被试目录
for sub_idx = 1:length(subject_names)
    current_subject = subject_names{sub_idx}; % 当前被试目录名，如'S001'
    subject_input_path = fullfile(input_root, current_subject);
    subject_output_path = fullfile(output_root, current_subject);
    
    % 创建输出目录（如果不存在）
    if ~exist(subject_output_path, 'dir')
        mkdir(subject_output_path);
        fprintf('创建输出目录: %s\n', subject_output_path);
    end
    
    % 获取该被试目录下的所有EDF文件
    edf_files = dir(fullfile(subject_input_path, '*.edf'));
    fprintf('正在处理被试 %s，共有 %d 个EDF文件\n', current_subject, length(edf_files));
    
    %% 4. 遍历处理当前被试下的每个EDF文件
    for file_idx =  1:length(edf_files)
        edf_file_name = edf_files(file_idx).name;
        edf_file_path = fullfile(subject_input_path, edf_file_name);
        
        fprintf('开始处理文件: %s\n', edf_file_name);
        
        try
            % 4.1 读取EDF文件
            eeg_struct = pop_biosig(edf_file_path);
            
            % 4.2 提取数据和相关参数
            fs = eeg_struct.srate; % 获取原始采样率
            
            nWant = 64;
            nChans = size(eeg_struct.data,1);
            EEG_channels = 1:min(nWant, nChans);
            EEG_data = eeg_struct.data(EEG_channels, :); % 提取指定通道的数据
            
            % 4.3 使用EEMD进行伪迹去除
            [~, channels_to_denoise] = ismember(target_channels, PHYSIONETMI_CHANNEL_LIST);
            channels_to_denoise = channels_to_denoise(channels_to_denoise > 0);
            
            nTotalChans = size(EEG_data, 1);
            channels_to_denoise = channels_to_denoise(channels_to_denoise <= nTotalChans);
            
            fprintf('选定 %d 个高风险通道进行EEMD去噪：\n', length(channels_to_denoise));
            for i = 1:length(channels_to_denoise)
                fprintf('%s ', PHYSIONETMI_CHANNEL_LIST{channels_to_denoise(i)});
            end
            fprintf('\n');
            
            % 4.4 仅对高风险通道进行EEMD去噪
            EEG_data_cleaned = EEG_data; % 初始化输出数据
            
            if ~isempty(channels_to_denoise)
                fprintf('开始EEMD去噪...\n');
                % 提取需要去噪的通道数据
                noisy_data = EEG_data(channels_to_denoise, :);
                % 应用EEMD去噪
                denoised_data = multichannel_eemd_denoise(noisy_data, fs, Nstd, NE, corr_threshold);
                % 将去噪后的数据放回原位
                EEG_data_cleaned(channels_to_denoise, :) = denoised_data;
            end
            
            clean_channels = setdiff(1:size(EEG_data, 1), channels_to_denoise);
            fprintf('其余 %d 个通道保留原始信号。\n', length(clean_channels));
            
%             EEG_data_cleaned = multichannel_eemd_denoise(EEG_data, fs, target_fs, Nstd, NE, corr_threshold);
            
            % 4.4 创建清理后的EEG结构体
            eeg_cleaned = eeg_struct; % 复制原结构体
            eeg_cleaned.data = EEG_data_cleaned; % 替换为清理后的数据
           
            
            % 4.5 定义输出路径并保存
            output_edf_path = fullfile(subject_output_path, edf_file_name);
           pop_writeeeg(eeg_cleaned, output_edf_path, 'TYPE', 'EDF');
            
            fprintf('成功处理并保存: %s\n', output_edf_path);
            
        catch ME
            % 处理过程中发生错误，记录错误信息但继续处理其他文件
            fprintf('错误处理文件 %s: %s\n', edf_file_name, ME.message);
            continue; % 继续处理下一个文件
        end
    end
end

fprintf('所有文件处理完成！\n');


