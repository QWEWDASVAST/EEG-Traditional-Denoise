clear;clc;
%% 准备文件路径
data_path = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\TSUBenchmark\Raw';
new_save_path = 'E:\workspace\Pycharmworkspace\EEGPT-main\datasets\TSUBenchmark\processing\ICA_EO_BSP';

% 确保输出目录存在
if ~isfolder(new_save_path)
    mkdir(new_save_path);
end

%% 获取所有.mat文件列表
% 使用dir函数获取目录下所有.mat文件 
mat_files = dir(fullfile(data_path, '*.mat'));
fprintf('找到%d个MAT文件需要处理\n', length(mat_files));

fs = 250;

%% 加载电极位置模板
load('biosemi_template.mat');          
topoplot([], locs, 'electrodes','ptslabels','plotdisk','on');

EEG_channels = 1:64;

%% 循环处理每个文件
for file_idx = 1:length(mat_files)
    try
        % 获取当前文件名
        data_file = mat_files(file_idx).name;
        fid = fullfile(data_path, data_file);
        
        fprintf('\n正在处理文件 %d/%d: %s\n', file_idx, length(mat_files), data_file);
        
        % 加载数据
        load(fid);
        
        fs = 250;
        
        %% 参数初始化
        [nChan, nTime, nFreq, nBlock] = size(data);
        nEpoch = nFreq * nBlock;               % 总 epoch 数
        
        % 初始化去噪后的4维数据数组
        data_denoised = zeros(nChan, nTime, nFreq, nBlock);
        
        %% 按epoch逐个进行ICA去伪迹处理
        fprintf('  开始逐个epoch处理，总共%d个epoch...\n', nEpoch);
        
        % 循环处理每个epoch
        epoch_count = 0;
        for freq_idx = 1:nFreq
            for block_idx = 1:nBlock
                epoch_count = epoch_count + 1;
                fprintf('  正在处理epoch %d/%d (频率%d/%d, 区块%d/%d)...\n', ...
                    epoch_count, nEpoch, freq_idx, nFreq, block_idx, nBlock);
                
                % 提取当前epoch的数据 
                current_epoch = data(:, :, freq_idx, block_idx);
                
                % 对当前epoch进行ICA去伪迹
                EEG_clean_epoch = IC_label_artifact_removal(current_epoch, fs, 826, locs, EEG_channels);
                
                % 验证处理后的数据维度
                [clean_nChan, clean_nTime] = size(EEG_clean_epoch);
                if clean_nChan ~= nChan || clean_nTime ~= nTime
                    fprintf('    警告: 处理后的epoch维度变化 (%d×%d -> %d×%d)，进行自动调整\n', ...
                        nChan, nTime, clean_nChan, clean_nTime);
                    % 如果维度不匹配，尝试调整或使用原始数据
                    if clean_nChan == nChan && clean_nTime <= nTime
                        % 如果通道数相同但时间点减少，用零填充
                        temp_data = zeros(nChan, nTime);
                        temp_data(:, 1:min(nTime, clean_nTime)) = EEG_clean_epoch(:, 1:min(nTime, clean_nTime));
                        EEG_clean_epoch = temp_data;
                    else
                        fprintf('    错误: 无法调整维度，使用原始epoch数据\n');
                        EEG_clean_epoch = current_epoch;
                    end
                end
                
                % 将处理后的epoch数据存回4维数组
                data_denoised(:, :, freq_idx, block_idx) = EEG_clean_epoch;
            end
        end
        
        fprintf('  所有epoch处理完成！\n');
        
        %% 7. 保存处理后的数据
        full_file_path = fullfile(new_save_path, data_file);
        % 加载时获取所有变量
        file_vars = load(fid);
        vars_list = fieldnames(file_vars);
        
        % 更新需要修改的变量
        file_vars.data = data_denoised;  % 替换为去噪后的数据
        
        % 保存所有变量（包括原始文件中可能存在的其他变量）
        save(full_file_path, '-struct', 'file_vars');
        fprintf('   已保留原始文件的所有变量结构，并更新了数据\n');
        
    catch ME
        fprintf('错误: 处理文件 %s 时发生错误: %s\n', data_file, ME.message);
        fprintf('错误堆栈:\n');
        for k = 1:length(ME.stack)
            fprintf('  %s (%d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        continue; % 跳过当前文件，继续处理下一个
    end
end

fprintf('\n批量处理完成！已处理 %d/%d 个文件\n', file_idx, length(mat_files));