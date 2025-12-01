clear; clc;

%%自动批量处理设置
input_folder = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\SEED\Preprocessed_EEG';
output_folder = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\SEED\EEMD';

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

% 准备电极位置
locs = readlocs('channel_62_pos.locs');
if size(locs, 1) > 62
    locs = locs(1:62, :);
    fprintf('电极位置已调整为62通道\n');
end

% EEMD参数设置 
Nstd = 0.2;          % 添加噪声的标准差比率
NE = 50;             % EEMD集成次数
corr_threshold = 0.1; % IMF相关系数阈值

% 为62通道SEED数据定义需要去噪的高风险通道
% 根据标准62通道系统定义高风险通道（前额和颞叶区域）
target_channels = {'Fp1', 'Fpz', 'Fp2', 'AF3', 'AF4', 'F7', 'F8', 'T7', 'T8'}; % 核心8通道

% 获取62通道的所有电极名称（从locs中提取）
if isstruct(locs)
    all_channel_names = {locs.labels};
else
    all_channel_names = locs;
end

fprintf('电极系统包含 %d 个通道\n', length(all_channel_names));

%% 批量处理所有文件
for file_idx = 1:length(file_list)
    current_file = file_list{file_idx};
    fprintf('\n正在处理文件 %d/%d: %s\n', file_idx, length(file_list), current_file);
    
    % 加载当前数据文件
    input_path = fullfile(input_folder, current_file);
    loaded_data = load(input_path);
    
    % 自动检测变量名模式
    field_names = fieldnames(loaded_data);
    fprintf('  文件中的字段: %s\n', strjoin(field_names, ', '));
    
    % 使用正则表达式检测变量名模式)
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
    
    % EEMD
    eeg_epochs_cell = cell(1, num_datasets);
    for i = 1:num_datasets
        fprintf('  处理数据集 %d/%d...\n', i, num_datasets);
        current_eeg_data = eeg_data_cell{i};
        
        % 获取目标通道的索引 - 与您提供的代码完全一致
        [~, channels_to_denoise] = ismember(target_channels, all_channel_names);
        channels_to_denoise = channels_to_denoise(channels_to_denoise > 0);

        nTotalChans = size(current_eeg_data, 1);
        channels_to_denoise = channels_to_denoise(channels_to_denoise <= nTotalChans);

        fprintf('    选定 %d 个高风险通道进行EEMD去噪：', length(channels_to_denoise));
        for j = 1:length(channels_to_denoise)
            fprintf('%s ', all_channel_names{channels_to_denoise(j)});
        end
        fprintf('\n');
        
        % 提取需要去噪的通道数据
        noisy_data = current_eeg_data(channels_to_denoise, :);
        
        % 应用EEMD去噪 - 直接调用您已有的multichannel_eemd_denoise函数
        denoised_data = multichannel_eemd_denoise(noisy_data, 200, Nstd, NE, corr_threshold);
        
        % 将去噪后的数据放回原位
        current_eeg_data(channels_to_denoise, :) = denoised_data;
    
        eeg_epochs_cell{i} = current_eeg_data;
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