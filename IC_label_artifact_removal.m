function [ EEG_data_processed] = IC_label_artifact_removal( EEG_data, fs, order,locs,EEG_channels)
addpath(genpath('E:\workspace\Pycharmworkspace\EEGdenoiseNet-master\code\preprocessing\lib'))
addpath(genpath('E:\workspace\Pycharmworkspace\EEGdenoiseNet-master\code\preprocessing\EEG\fastica'))

%% 2. rereference
low_cutoff = 1;    
% 根据采样率自动设定高通截止频率80；79,63； 
if fs == 160
    high_cutoff = 79;      
elseif fs == 128
    high_cutoff = 63;     
else
    high_cutoff = 80;     
end
EEG_data = filter_data(EEG_data, fs, low_cutoff, high_cutoff, order);% 通滤波（1-80 Hz），滤除极低频漂移和高频噪声（如肌电）

%% 4. ICA decomposition using fast ICA
% 3.1 pca withening
[dewhitening, score, latent]=pca(EEG_data','Centered','off');
whitening=inv(dewhitening);
PC=whitening*EEG_data;
% 计算主成分的方差贡献率，并根据阈值（0.01%的平均方差）确定保留的主成分数量 (n_pcs)。这用于降维，减少计算量并避免过拟合。
latent=100*latent/sum(latent);
n_pcs=sum(latent>0.01*mean(latent));
retained_variance=sum(latent(1:n_pcs));
[IC, mixing, unmixing]=fastica(PC(1:n_pcs, :),'approach','defl','g','tanh','maxNumIterations',500); 
IC=unmixing*PC(1:n_pcs,:);
mixing_matrix = dewhitening(:,1:n_pcs)*mixing;
unmixing_matrix = unmixing*whitening(1:n_pcs,:);

%% 4. Use ICLabel
% prepare EEGlab struct （为ICLabel创建EEGLAB标准数据结构）
data_len = size(EEG_data,2);
time_len = data_len/fs;
time_axis = linspace(0,time_len, data_len);
EEG_struct.times = time_axis;
EEG_struct.data = EEG_data;
EEG_struct.chanlocs = locs;
EEG_struct.srate = fs;
EEG_struct.trials = 1;
EEG_struct.pnts = data_len;
EEG_struct.icawinv = mixing_matrix;
EEG_struct.icaweights = unmixing_matrix;
EEG_struct.icaact = unmixing_matrix*EEG_data;
EEG_struct.icachansind = EEG_channels;
EEG_struct.ref = 'averef';
% classification 
EEG_struct= iclabel(EEG_struct); % 调用ICLabel进行分类
class_prob = EEG_struct.etc.ic_classification.ICLabel.classifications;
class_type = EEG_struct.etc.ic_classification.ICLabel.classes;
%% 7. reconstruct cleaned EEG and save to D
ic_num = size(class_prob,1);
for iter_ic = 1:ic_num
        [max_prob_val(1, iter_ic), class_index(1, iter_ic)] = max(class_prob(iter_ic, :));
end
class_label.brain_ic_index = sort(find(class_index == 1));
class_label.muscle_ic_index = sort(find(class_index == 2));
class_label.eye_ic_index = sort(find(class_index == 3));
class_label.heart_ic_index = sort(find(class_index == 4));
class_label.line_noise_ic_index = sort(find(class_index == 5));
class_label.channel_noise_ic_index = sort(find(class_index == 6));
class_label.other_noise_ic_index = sort(find(class_index == 7)); % what is other, should I remove it?
brain_ic_num = length(class_label.brain_ic_index);
if(brain_ic_num == 0)
        warning('No brain IC recovered! Please check')
end
    
bad_ic_index = unique(sort([class_label.muscle_ic_index, class_label.eye_ic_index, class_label.heart_ic_index, class_label.line_noise_ic_index, class_label.channel_noise_ic_index]));%, class_label.other_noise_ic_index
    
EEG_data_processed = EEG_data - mixing_matrix(:, bad_ic_index)*IC(bad_ic_index,:);

% %将重建后的信号转换为双精度格式，并重采样到目标采样率 ，以统一采样率便于后续分析。
% EEG_data_processed  = double(EEG_data );
% 
% epoch_length_samples = fs; % 1秒
% t_num = size(EEG_data_processed,2);
% % good_data_num = idivide(t_num,int32(1000));
% good_data_num = floor(t_num / epoch_length_samples);
% 
% 
% %分段与Epoch筛选
% %将连续的EEG信号分割成1秒长的epochs。
% EEG_epochs = [];
% 
% for iter_cut = 1:good_data_num
%     start = (iter_cut - 1) * epoch_length_samples + 1;
%     end_idx = start + epoch_length_samples - 1;
%     EEG_epochs = [EEG_epochs; EEG_data_processed(:, start:end_idx)];
% end
% 
% %% 4. remove bad epochs according to PSD template
% %基于功率谱密度 (PSD) 模板匹配进一步筛选epoch。计算每个epoch的PSD与一个预定义的模板PSD
% %（通常代表“良好”EEG的频谱特征，如低频能量较高）的相关系数，剔除相关系数过低（<0.8）的epoch。
% % frequencies = 1:120;
% % template_PSD = exp(-1*frequencies/20);
% epoch_num = size(EEG_epochs, 1);
% corr_threshold = 0.65;
% frequencies = 1:min(120,floor(fs/2));
% template_PSD = exp(-frequencies/20);
% template_PSDs = repmat(template_PSD', 1, epoch_num);
% 
% [PSDs, f] = pwelch(EEG_epochs', [], [], frequencies, fs);
% 
% PSD_corr = corr(PSDs, template_PSDs);
% PSD_corr = PSD_corr(:,1);
% EEG_epochs(PSD_corr < corr_threshold, :) = []; % remove those are not look like template
% PSDs(:, PSD_corr < corr_threshold) = [];
% 
% %计算每个epoch的高频（40-120 Hz）功率占总功率的比例。
% %比例过高（>0.2）可能表明该epoch仍残留肌电伪迹（肌电信号具有丰富的高频成分），因此将其剔除。
% 
% high_freq = 40:min(120, floor(fs/2));
% power_all = sum(PSDs, 1);
% high_power = sum(PSDs(high_freq,:), 1);
% ratio = high_power./power_all;
% EEG_epochs(ratio > 0.2,: ) = []; % when there's too much high frequency, it might be contaminated by EMG
% PSDs(:, ratio >0.2) = [];
% 
% %% 4. remove bad EMGs according to PSD
% % ratio_threshold = 50; % low_freq/high_freq < 1, the more higher frequency, the better the EMG signal
% % % frequencies = 1:200;
% % % low_freq = 1:30;
% % % high_freq = 40:200;
% % frequencies = 1:floor(fs/2);
% % low_freq = 1:min(30, floor(fs/2));
% % high_freq = 40:min(200, floor(fs/2));
% % 
% 
% % [PSDs, f] = pwelch(EEG_epochs', [], [], frequencies, fs);
% % low_powers = mean(PSDs(low_freq,:), 1);
% % high_powers = mean(PSDs(high_freq,:), 1);
% % ratios = low_powers./high_powers;
% % %计算每个epoch的低频（1-30 Hz）与高频（40-200 Hz）平均功率的比值。
% % %比值过低（<70）可能表示高频噪声（如肌电）污染严重，因此将其剔除。
% % EEG_epochs(ratios<ratio_threshold, :) = []; % remove those are not look like template
% 
% %% 修改后的Epoch重组部分 - 使用NaN填充被删除的epoch
% 
% % 获取原始数据的通道数和时间点数
% [n_channels_original, n_timepoints_original] = size(EEG_data);
% 
% % 检查EEG_epochs是否为空
% if isempty(EEG_epochs)
%     warning('所有epoch都被剔除，返回原始ICA处理后的数据');
%     return;
% end
% 
% % 获取筛选后epoch的信息
% [epoch_rows, epoch_samples] = size(EEG_epochs);
% 
% % 计算完整epoch的数量（考虑可能的不完整情况）
% n_complete_epochs = floor(epoch_rows / n_channels_original);
% 
% if n_complete_epochs == 0
%     warning('没有完整的epoch保留，返回原始数据');
%     return;
% end
% 
% % 如果有不完整的epoch，截取完整部分
% if mod(epoch_rows, n_channels_original) ~= 0
%     epoch_rows = n_complete_epochs * n_channels_original;
%     EEG_epochs = EEG_epochs(1:epoch_rows, :);
%     fprintf('截断不完整的epoch，保留%d个完整epoch\n', n_complete_epochs);
% end
% 
% %% 方案：使用NaN填充被删除的epoch，保持原始数据形状
% % 创建与原始数据相同大小的NaN矩阵
% EEG_data_reconstructed = NaN(n_channels_original, n_timepoints_original);
% 
% % 计算原始数据中总共可以分成多少个完整的epoch
% total_possible_epochs = floor(n_timepoints_original / epoch_length_samples);
% 
% % 我们需要确定哪些epoch被保留了，哪些被删除了
% % 由于删除操作是在循环中进行的，我们需要重构索引映射
% 
% % 方法：按顺序将保留的epoch放置到原始位置，被删除的位置保持NaN
% current_epoch_in_original = 1;
% current_epoch_in_remaining = 1;
% 
% % 记录每个原始epoch位置的保留状态
% epoch_preservation_status = false(1, total_possible_epochs);
% 
% % 由于我们不知道具体哪些epoch被删除了，采用顺序映射的方法
% % 假设保留的epoch按时间顺序排列
% for epoch_idx = 1:total_possible_epochs
%     % 计算当前epoch在原始数据中的位置
%     start_sample = (epoch_idx - 1) * epoch_length_samples + 1;
%     end_sample = min(epoch_idx * epoch_length_samples, n_timepoints_original);
%     
%     % 检查是否还有保留的epoch需要放置
%     if current_epoch_in_remaining <= n_complete_epochs && ...
%        start_sample <= n_timepoints_original && ...
%        end_sample <= n_timepoints_original
%         
%         % 计算当前保留的epoch在EEG_epochs中的位置
%         start_row = (current_epoch_in_remaining - 1) * n_channels_original + 1;
%         end_row = current_epoch_in_remaining * n_channels_original;
%         
%         % 提取当前保留的epoch数据
%         current_epoch_data = EEG_epochs(start_row:end_row, :);
%         
%         % 将保留的epoch数据放入对应位置
%         EEG_data_reconstructed(:, start_sample:end_sample) = current_epoch_data;
%         
%         % 标记这个epoch被保留了
%         epoch_preservation_status(epoch_idx) = true;
%         current_epoch_in_remaining = current_epoch_in_remaining + 1;
%     end
%     % 如果条件不满足，该位置保持NaN（表示epoch被删除）
% end
% 
% % 处理可能剩余的数据点（如果原始数据长度不是epoch长度的整数倍）
% remaining_samples = n_timepoints_original - total_possible_epochs * epoch_length_samples;
% if remaining_samples > 0 && current_epoch_in_remaining <= n_complete_epochs
%     % 如果还有保留的epoch且有余数样本，尝试放置
%     start_sample = total_possible_epochs * epoch_length_samples + 1;
%     end_sample = n_timepoints_original;
%     
%     start_row = (current_epoch_in_remaining - 1) * n_channels_original + 1;
%     end_row = current_epoch_in_remaining * n_channels_original;
%     
%     % 只放置实际存在的样本
%     available_samples = min(epoch_length_samples, remaining_samples);
%     if size(EEG_epochs, 2) >= available_samples
%         current_epoch_data = EEG_epochs(start_row:end_row, 1:available_samples);
%         EEG_data_reconstructed(:, start_sample:end_sample) = current_epoch_data;
%     end
% end
% 
% % 用重构的数据替换原来的EEG_data_processed
% EEG_data_processed = EEG_data_reconstructed;


    