clear; clc;

%% 1. 自动批量处理设置
input_folder = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\SEED\Preprocessed_EEG';
output_folder = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\SEED\WT';

% 创建输出文件夹
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    fprintf('创建输出目录: %s\n', output_folder);
end

% 获取需要处理的所有.mat文件列表（排除label.mat）
all_files = dir(fullfile(input_folder, '*.mat'));
file_list = {};
for i = 1:length(all_files)
    if ~strcmp(all_files(i).name, 'label.mat')
        file_list{end+1} = all_files(i).name;
    end
end

fprintf('找到 %d 个需要处理的文件\n', length(file_list));

% 设置EEG通道
EEG_channels = 1:62;

% 小波阈值去噪参数
wname = 'db4';        % 小波名称，db4小波在EEG处理中常用
level = 3;            % 小波分解水平

% 准备电极位置
locs = readlocs('channel_62_pos.locs');
if size(locs, 1) > 62
    locs = locs(1:62, :);
    fprintf('电极位置已调整为62通道\n');
end

%%批量处理所有文件
for file_idx = 1:length(file_list)
    current_file = file_list{file_idx};
    fprintf('\n正在处理文件 %d/%d: %s\n', file_idx, length(file_list), current_file);
    
    % 加载当前数据文件
    input_path = fullfile(input_folder, current_file);
    loaded_data = load(input_path);
    
    % 自动检测变量名模式
    field_names = fieldnames(loaded_data);
    fprintf('  文件中的字段: %s\n', strjoin(field_names, ', '));
    
    % 使用正则表达式检测变量名模式
    eeg_pattern = '^([a-zA-Z]+)_eeg\d+$';
    eeg_fields = {};
    
    for i = 1:length(field_names)
        if ~isempty(regexp(field_names{i}, eeg_pattern, 'once'))
            eeg_fields{end+1} = field_names{i};
        end
    end
    
    if isempty(eeg_fields)
        error('在文件 %s 中未找到符合*_eeg*模式的变量', current_file);
    end
    
    % 提取前缀并排序字段
    prefix = regexp(eeg_fields{1}, '^([a-zA-Z]+)_eeg', 'tokens');
    prefix = prefix{1}{1};
    
    % 按数字排序eeg字段
    eeg_numbers = cellfun(@(x) str2double(regexp(x, '\d+$', 'match')), eeg_fields);
    [~, sort_idx] = sort(eeg_numbers);
    eeg_fields_sorted = eeg_fields(sort_idx);
    
    num_datasets = length(eeg_fields_sorted);
    fprintf('  检测到变量模式: %s_eeg*, 共%d个数据集\n', prefix, num_datasets);
    
    % 提取数据到元胞数组中
    eeg_data_cell = cell(1, num_datasets);
    for i = 1:num_datasets
        var_name = eeg_fields_sorted{i};
        eeg_data_cell{i} = loaded_data.(var_name);
        fprintf('  加载 %s: 尺寸 %dx%d\n', var_name, size(eeg_data_cell{i}, 1), size(eeg_data_cell{i}, 2));
    end
    
    % 使用ICA和IClabel方法处理所有数据集
    eeg_epochs_cell = cell(1, num_datasets);
    for i = 1:num_datasets
        fprintf('  处理数据集 %d/%d...\n', i, num_datasets);
        current_eeg_data = eeg_data_cell{i};
        eeg_epochs_cell{i} = WT_artifact_removal(current_eeg_data, EEG_channels, wname, level);
    end
    
    % 保存处理后的数据，保持与原文件一致的结构和文件名
    output_file = fullfile(output_folder, current_file);

    % 创建结构体来存储所有处理后的数据
    output_data = struct();
    for i = 1:num_datasets
        var_name = eeg_fields_sorted{i};
        output_data.(var_name) = eeg_epochs_cell{i};
    end

    % 使用save函数保存，避免v7.3格式
    save(output_file, '-struct', 'output_data');

    
    fprintf('文件已成功保存: %s\n', output_file);
    fprintf('保存的变量: %s\n', strjoin(eeg_fields_sorted, ', '));
end

fprintf('\n批量处理完成! 共处理 %d 个文件\n', length(file_list));


%% 小波阈值去噪函数
function EEG_data_cleaned = WT_artifact_removal(EEG_data, EEG_channels, wname, level)
    % 使用小波阈值进行伪迹去除（使用标准小波阈值去噪方法）
    fprintf('正在进行小波阈值去噪...\n');
    
    % 提取指定通道的数据
    EEG_data_selected = EEG_data(EEG_channels, :);
    EEG_data_cleaned = zeros(size(EEG_data_selected));
    
    % 对每个通道分别进行小波阈值去噪
    for ch = 1:size(EEG_data_selected, 1)
        channel_data = EEG_data_selected(ch, :);
        N = length(channel_data);
        
        % 使用标准小波阈值去噪方法
        % 小波分解
        [c, l] = wavedec(channel_data, level, wname);
        
        % 对每个细节系数进行软阈值处理
        for i = 1:level
            start_index = sum(l(1:i)) + 1;
            end_index = sum(l(1:i+1));
            
            % 使用通用阈值公式
            thr = 0.5 * sqrt(2 * log(N)); % 阈值0.5*sqrt(2*log(N))
            
            % 应用软阈值处理
            c(start_index:end_index) = wthresh(c(start_index:end_index), 's', thr); 
        end
        
        % 将去噪后的信号进行小波重构
        cleaned_channel = waverec(c, l, wname);
        
        EEG_data_cleaned(ch, :) = cleaned_channel;
    end
end