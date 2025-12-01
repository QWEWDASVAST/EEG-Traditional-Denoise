clear; clc;
addpath(genpath('E:\workspace\Pycharmworkspace\EEGdenoiseNet-master\code\preprocessing\EEG\waveletTh_wp'))

%% 1. 定义路径和参数
input_root = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\PhysioNetMI\files\';
output_root = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\PhysioNetMI\processing\WT\pycache_data\';

% 加载电极位置模板（假设所有文件共用同一套电极位置）
load('biosemi_template.mat'); % 确保此文件在MATLAB路径中

% 定义要处理的通道
EEG_channels = 1:64;
% 小波阈值去噪参数
wname = 'db4';        % 小波名称，db4小波在EEG处理中常用
level = 3;            % 小波分解水平


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
            
            % 4.3 使用小波阈值进行伪迹去除（使用标准小波阈值去噪方法）
            fprintf('正在进行小波阈值去噪...\n');
            EEG_data_cleaned = zeros(size(EEG_data));
            
            % 对每个通道分别进行小波阈值去噪[6,7](@ref)
            for ch = 1:size(EEG_data, 1)
                channel_data = EEG_data(ch, :);
                N = length(channel_data);
                
                % 使用标准小波阈值去噪方法[6,7](@ref)
                % 小波分解
                [c, l] = wavedec(channel_data, level, wname);
                
                % 对每个细节系数进行软阈值处理
                for i = 1:level
                    start_index = l(i) + 1;
                    end_index = l(i + 1);
                    
                    % 使用通用阈值公式
                    thr = 0.5 * sqrt(2 * log(N)); % 阈值0.5*sqrt(2*log(N))
                    
                    % 应用软阈值处理
                    c(start_index:end_index) = wthresh(c(start_index:end_index), 's', thr); 
                end
                
                % 将去噪后的信号进行小波重构
                cleaned_channel = waverec(c, l, wname);
                
                EEG_data_cleaned(ch, :) = cleaned_channel;
            end
            
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


