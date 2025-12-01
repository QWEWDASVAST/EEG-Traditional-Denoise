clear; clc;

%% 1. 定义路径和参数
input_root = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\PhysioNetMI\files\';
output_root = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\PhysioNetMI\processing\ICA_EO_BSP1\pycache_data\';

% 加载电极位置模板
load('biosemi_template.mat');

% 定义要处理的通道
EEG_channels = 1:64;

%% 2. 获取所有被试目录
% 使用dir获取所有以'S'开头的文件夹
subject_folders = dir(fullfile(input_root, 'S*'));
subject_folders = subject_folders([subject_folders.isdir]); 
subject_names = {subject_folders.name}; 
fprintf('找到 %d 个被试目录\n', length(subject_names));

%% 3. 遍历每个被试目录
for sub_idx = 1:length(subject_names)
    current_subject = subject_names{sub_idx}; 
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
    for file_idx = 1:length(edf_files)
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
            
            % 4.3 使用ICA和ICLabel进行伪迹去除
            EEG_data_cleaned = IC_label_artifact_removal(EEG_data, fs, 528, locs, EEG_channels);
            
            % 4.4 创建清理后的EEG结构体
            eeg_cleaned = eeg_struct;
            eeg_cleaned.data = EEG_data_cleaned;
                     
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


