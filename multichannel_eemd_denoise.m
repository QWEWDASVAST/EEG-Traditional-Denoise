function [denoised_data] = multichannel_eemd_denoise(EEG_data, fs, Nstd, NE, corr_threshold)
    % 函数：对多通道EEG数据进行EEMD去噪
    % 输入：
    %   EEG_data: 多通道EEG数据 (channels x samples)
    %   fs: 原始采样率
    %   Nstd: 添加噪声的标准差比率
    %   NE: EEMD集成次数
    %   corr_threshold: IMF相关系数阈值
    % 输出：
    %   denoised_data: 去噪后的EEG数据
    
    [num_channels, num_samples] = size(EEG_data);
    denoised_data = zeros(num_channels, num_samples);
    
    fprintf('开始多通道EEMD去噪，使用相关系数阈值: %.3f\n', corr_threshold);
    fprintf('参数设置 - Nstd: %.2f, NE: %d\n', Nstd, NE);
    
    % 检查并行池状态，如未启动则启动并行池
    if isempty(gcp('nocreate'))
        try
            parpool('local'); % 使用默认worker数
            fprintf('已启动并行池\n');
        catch ME
            fprintf('无法启动并行池，将使用串行处理: %s\n', ME.message);
        end
    else
        pool = gcp;
        fprintf('并行池已存在，%d个worker可用\n', pool.NumWorkers);
    end
    
    % 使用parfor并行处理每个通道
    parfor ch = 1:num_channels
        channel_start_time = tic;
        
        try
            % 提取当前通道数据
            channel_data = EEG_data(ch, :);
            
            % 执行EEMD分解
            allmode = eemd(channel_data, Nstd, NE);
            
            % 获取IMF数量
            num_components = size(allmode, 2);
            num_imfs = num_components - 2; % 减去原始信号和残差
            
            if num_imfs <= 0
                warning('通道 %d: 无法分解出有效的IMF，使用原始信号。', ch);
                denoised_channel = channel_data;
            else
                % 计算每个IMF与原始信号的相关系数
                correlations = zeros(1, num_imfs);
                original_signal = allmode(:, 1);
                
                for i = 1:num_imfs
                    imf_signal = allmode(:, i+1);
                    corr_matrix = corrcoef(original_signal, imf_signal);
                    correlations(i) = abs(corr_matrix(1, 2));
                end
                
                % 根据相关系数阈值选择要保留的IMF
                imfs_to_keep = find(correlations > corr_threshold);
                
                if isempty(imfs_to_keep)
                    % 如果所有IMF相关系数都太低，保留相关性最高的前几个
                    [~, sorted_indices] = sort(correlations, 'descend');
                    imfs_to_keep = sorted_indices(1:min(3, num_imfs));
                end
                
                % 重构去噪后的信号
                imf_indices = imfs_to_keep + 1;
                residual_index = size(allmode, 2);
                denoised_channel = sum(allmode(:, [imf_indices, residual_index]), 2)';
                
                % 记录处理信息（parfor中fprintf是安全的）
                fprintf('通道 %d: 保留 %d/%d 个IMF，最高相关系数: %.3f\n', ...
                    ch, length(imfs_to_keep), num_imfs, max(correlations));
            end
            
            % 存储去噪后的通道数据
            denoised_data(ch, :) = denoised_channel;
            
            channel_time = toc(channel_start_time);
            fprintf('通道 %d 处理完成，耗时: %.2f秒\n', ch, channel_time);
            
        catch ME
            % 如果处理失败，使用原始信号并记录错误
            fprintf('通道 %d 处理失败: %s，使用原始信号\n', ch, ME.message);
            denoised_data(ch, :) = EEG_data(ch, :);
        end
    end
    
    fprintf('多通道EEMD去噪完成\n');
end