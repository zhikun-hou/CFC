
% 小波熵
% 由于小波熵基于离散小波变换，而非连续小波变换，因此不能直接使用ft_freqanalysis的wavelet来计算
% Wavelet entropy: a new tool for analysis of short duration brain electrical signals

function [Entropy] = ft_entropyWT(cfg,FT)


    % 预处理 ==============================================================

    FT  = ft_checkdata(FT,'datatype','raw','feedback','yes');
    
    cfg.trials  = ft_getopt(cfg,'trials',1:numel(FT.trial));
    cfg.channel = ft_getopt(cfg,'channel',1:numel(FT.label));
    cfg.latency = 'minperiod'; % 使用所有试次中最短的试次对所有试次进行截取，避免长短不一致的问题

    % 采样率通常在1024左右，那么奈奎斯特频率为512Hz
    % 7阶时划分了0-4Hz和4-8Hz，通常够用了
    % 0-4Hz 慢波  4-8Hz theta波  8-16Hz alpha波  16-32Hz beta波
    % 32-64Hz LowGamma  64-128Hz  128-256Hz  256-512Hz
    % 注意，七阶划分出了八个频段
    cfg.order = ft_getopt(cfg,'order',7);
    
    FT = ft_selectdata(cfg,FT);

    % 检查公共设置 =========================================================

    cfg.wavelet = ft_getopt(cfg,'wavelet','db4'); % 所使用的wavelet的name，见wfilters文档的wname

    cfg.toi = ft_getopt(cfg,'toi','all');
    if(~strcmp(cfg.toi,'all') && ~isnumeric(cfg.toi))
        ft_error("toi应当为all或者(1,N)或者(N,1)，N>=2");
    end
    
    % 如果toi的N>2，就需要指定twin，注意是半长而不是全长
    cfg.twin = ft_getopt(cfg,'twin',0.1);

    % 注意，toi指定的是时间窗的中心时刻
    if(isnumeric(cfg.toi) && ~issorted(cfg.toi))
        ft_error("手动指定toi时，应当单调递增");
    end

    cfg.baseline = ft_getopt(cfg,'baseline','none');
    if(isnumeric(cfg.baseline))
        if(numel(cfg.baseline)~=2)
            ft_error("baseline应当为(1,2)或者(2,1)");
        end
        if(~issorted(cfg.baseline))
            ft_error("手动指定baseline时，应当单调递增");
        end
    end
    cfg.baselinetype = ft_getopt(cfg,'baselinetype','absolute');


    % 如果没有指定是否可视化，并且nargout为0，那就自动执行可视化
    cfg.visualize = ft_getopt(cfg,'visualize',nargout==0);

    % 执行 ================================================================

    % 共有Total Wavelet Entropy、Relative Wavelet Entropy、
    % Time Evolution of WE/RWE、Event-Related WE 五种模式
    
    switch(cfg.baselinetype)
        case {'absolute','total','no','none'}
            if(strcmp(cfg.toi,'all') || numel(cfg.toi)==2)
                Entropy = getTotalWaveletEntropy(cfg,FT);
            else
                Entropy = getTimeEvolutionTotalWaveletEntropy(cfg,FT);
            end
        case {'relative'}
            if(strcmp(cfg.toi,'all') || numel(cfg.toi)==2)
                Entropy = getRelativeWaveletEntropy(cfg,FT);
            else
                Entropy = getTimeEvolutionRelativeWaveletEntropy(cfg,FT);
            end
        case {'relchange'}
            Entropy     = getEventRelatedTotalWaveletEntropy(cfg,FT);
        otherwise
            ft_error("未知的baselinetype");
    end


end




function [E] = getEventRelatedTotalWaveletEntropy(cfg,FT)

    N_trial   = numel(FT.trial);
    N_channel = numel(FT.label);
    N_order = cfg.order;
    N_window = numel(cfg.toi);

    % 离散小波变换
    % 注意，DWT是霍夫曼二叉树，从第n层的approx分出n+1层的approx和detail
    % 因此总能量=各层的detail+最后一层的approx
    Energy = cell(1,N_trial);
    Entropy = zeros(N_window,N_channel,N_trial);
    % 此外，需要注意：分解是从高频开始，低通/高通的分界线每次/=2
    % 第一次分解是Fs/4~Fs/2为高频、作为detail，0~Fs/4为低频
    ft_progress('init','etf');
    for i=1:N_trial
        N_time = size(FT.trial{i},2);

        ft_progress(i/N_trial,'正在处理试次(%d/%d)', i, N_trial);
            
        % 计算基线对应的索引
        Idx_baseline = FT.time{i}>=cfg.baseline(1) & FT.time{i}<=cfg.baseline(2);

        for j=1:N_channel
            Pow = zeros(N_time,N_order+1);
            Freq = zeros(1,N_order+2); % 因为每次/=2，所以这里存储频段的边界而非中心

            [C,L] = wavedec(FT.trial{i}(j,:),N_order,cfg.wavelet);
            
            % 分解时是从高频到低频，但我们存储时按照从低频到高频排列
            Freq(1) = FT.fsample / 2;
            for k=1:N_order
                idx = N_order+2-k; % 要存储的位置
                Pow(:,idx) = wrcoef('d',C,L,cfg.wavelet,k).^2; % 重建小波系数，然后计算功率

                Freq(k+1) = Freq(k) / 2;
            end
            Freq = flip(Freq); % 翻转，变成从低频到高频
            Pow(:,1) = wrcoef('a',C,L,cfg.wavelet,N_order).^2; % (times,bands)
            
            % 计算基线
            Pow_baseline     = sum(Pow(Idx_baseline,:),1);
            Energy{i,1}(:,j) = Pow_baseline; % (1,bands)=>(bands,channels)
            P_baseline       = Pow_baseline ./ sum(Pow_baseline);
            H_baseline       = -sum(P_baseline.*log(P_baseline));

            % 逐个计算时间窗
            for n=1:N_window
                Idx_win = FT.time{i}>=cfg.toi(n)-cfg.twin & FT.time{i}<=cfg.toi(n)+cfg.twin;
    
                % 计算能量
                Pow_toi      = sum(Pow(Idx_win,:),1);
                P_toi        = Pow_toi ./ sum(Pow_toi); % (1,bands)
                Energy{i,2}(n,:,j) = Pow_toi; % bands=orders+1
                
                % Relative Wavelet Energy=特定尺度下所有时刻的能量之和/总能量
                H_toi          = -sum(P_toi.*log(P_toi));
                Entropy(n,j,i) = (H_toi-H_baseline)/H_baseline;
            end
        end
        
    end
    ft_progress('close');


    % 输出 ================================================================

    E = [];
    E.method    = "EventRelatedTotalWaveletEntropy";
    E.channels  = cfg.channel;
    E.trials    = cfg.trials;
    E.times     = cfg.toi;
    E.freqs     = Freq;

    E.Energy  = Energy; % 每个时刻、每个尺度的能量
    E.Entropy = Entropy; 
    
    if(cfg.visualize)
        AVG = mean(E.Entropy,3);
        SE  = std(E.Entropy,[],3)/sqrt(N_trial);
        for j=1:N_channel
            patch('XData',[E.times flip(E.times)],'YData',[AVG(:,j)+SE(:,j);flip(AVG(:,j)-SE(:,j),1)],'FaceAlpha',0.1,'FaceColor','b','EdgeAlpha',0);
        end
        hold on;
        plot(E.times,AVG,'r');
    end

end



function [E] = getRelativeWaveletEntropy(cfg,FT)

    N_trial   = numel(FT.trial);
    N_channel = numel(FT.label);
    N_order = cfg.order;

    % 离散小波变换
    % 注意，DWT是霍夫曼二叉树，从第n层的approx分出n+1层的approx和detail
    % 因此总能量=各层的detail+最后一层的approx
    Energy = cell(N_trial,2);
    Entropy = zeros(N_channel,N_trial);
    % 此外，需要注意：分解是从高频开始，低通/高通的分界线每次/=2
    % 第一次分解是Fs/4~Fs/2为高频、作为detail，0~Fs/4为低频
    ft_progress('init','etf');
    for i=1:N_trial
        N_time = size(FT.trial{i},2);

        ft_progress(i/N_trial,'正在处理试次(%d/%d)', i, N_trial);

        % 计算时间对应的索引
        Idx_baseline = FT.time{i}>=cfg.baseline(1) & FT.time{i}<=cfg.baseline(2);
        if(strcmp(cfg.toi,'all'))
            Idx_toi  = true(1,N_time);
        else
            Idx_toi  = FT.time{i}>=cfg.toi(1) & FT.time{i}<=cfg.toi(2);
        end

        for j=1:N_channel
            Pow = zeros(N_time,N_order+1);
            Freq = zeros(1,N_order+2); % 因为每次/=2，所以这里存储频段的边界而非中心

            [C,L] = wavedec(FT.trial{i}(j,:),N_order,cfg.wavelet);
            
            % 分解时是从高频到低频，但我们存储时按照从低频到高频排列
            Freq(1) = FT.fsample / 2;
            for k=1:N_order
                idx = N_order+2-k; % 要存储的位置
                Pow(:,idx) = wrcoef('d',C,L,cfg.wavelet,k).^2; % 重建小波系数，然后计算功率

                Freq(k+1) = Freq(k) / 2;
            end
            Freq = flip(Freq); % 翻转，变成从低频到高频
            Pow(:,1) = wrcoef('a',C,L,cfg.wavelet,N_order).^2; % (times,bands)

            % 计算能量
            Pow_toi      = sum(Pow(Idx_toi,:),1);
            Pow_baseline = sum(Pow(Idx_baseline,:),1);

            Energy{i,1}(:,j) = Pow_baseline; % (1,bands)=>(bands,channels)
            Energy{i,2}(:,j) = Pow_toi; % bands=orders+1
            
            % Relative Wavelet Energy=特定尺度下所有时刻的能量之和/总能量
            P_toi      = Pow_toi ./ sum(Pow_toi); % (1,bands)
            P_baseline = Pow_baseline ./ sum(Pow_baseline);
            
            % Relative Wavelet Entropy
            Entropy(j,i) = sum(P_toi.*log(P_toi./P_baseline));

        end

    end
    ft_progress('close');

    % 输出 ================================================================

    E = [];
    E.method    = "RelativeWaveletEntropy";
    E.channels  = cfg.channel;
    E.trials    = cfg.trials;
    E.freqs     = Freq;

    E.Energy  = Energy; % 每个时刻、每个尺度的能量
    E.Entropy = Entropy; 
    
    if(cfg.visualize)
        figure("Name","RelativeWaveletEntropy");
        imagesc(E.trials,E.channels,E.Entropy);
        colorbar;
        clim([0 2]);
        colormap(jet);
    end

end



function [E] = getTimeEvolutionRelativeWaveletEntropy(cfg,FT)

    N_trial   = numel(FT.trial);
    N_channel = numel(FT.label);
    N_order = cfg.order;
    N_window = numel(cfg.toi);

    % 离散小波变换
    % 注意，DWT是霍夫曼二叉树，从第n层的approx分出n+1层的approx和detail
    % 因此总能量=各层的detail+最后一层的approx
    Energy = cell(1,N_trial);
    Entropy = zeros(N_window,N_channel,N_trial);
    % 此外，需要注意：分解是从高频开始，低通/高通的分界线每次/=2
    % 第一次分解是Fs/4~Fs/2为高频、作为detail，0~Fs/4为低频
    ft_progress('init','etf');
    for i=1:N_trial
        N_time = size(FT.trial{i},2);

        ft_progress(i/N_trial,'正在处理试次(%d/%d)', i, N_trial);
            
        % 计算基线对应的索引
        Idx_baseline = FT.time{i}>=cfg.baseline(1) & FT.time{i}<=cfg.baseline(2);

        for j=1:N_channel
            Pow = zeros(N_time,N_order+1);
            Freq = zeros(1,N_order+2); % 因为每次/=2，所以这里存储频段的边界而非中心

            [C,L] = wavedec(FT.trial{i}(j,:),N_order,cfg.wavelet);
            
            % 分解时是从高频到低频，但我们存储时按照从低频到高频排列
            Freq(1) = FT.fsample / 2;
            for k=1:N_order
                idx = N_order+2-k; % 要存储的位置
                Pow(:,idx) = wrcoef('d',C,L,cfg.wavelet,k).^2; % 重建小波系数，然后计算功率

                Freq(k+1) = Freq(k) / 2;
            end
            Freq = flip(Freq); % 翻转，变成从低频到高频
            Pow(:,1) = wrcoef('a',C,L,cfg.wavelet,N_order).^2; % (times,bands)
            
            % 计算基线
            Pow_baseline     = sum(Pow(Idx_baseline,:),1);
            Energy{i,1}(:,j) = Pow_baseline; % (1,bands)=>(bands,channels)
            P_baseline       = Pow_baseline ./ sum(Pow_baseline);

            % 逐个计算时间窗
            for n=1:N_window
                Idx_win = FT.time{i}>=cfg.toi(n)-cfg.twin & FT.time{i}<=cfg.toi(n)+cfg.twin;
    
                % 计算能量
                Pow_toi      = sum(Pow(Idx_win,:),1);
                P_toi        = Pow_toi ./ sum(Pow_toi); % (1,bands)
                Energy{i,2}(n,:,j) = Pow_toi; % bands=orders+1
                
                % Relative Wavelet Energy=特定尺度下所有时刻的能量之和/总能量
                Entropy(n,j,i) = sum(P_toi.*log(P_toi./P_baseline));
            end
        end
        
    end
    ft_progress('close');


    % 输出 ================================================================

    E = [];
    E.method    = "TimeEvolutionRelativeWaveletEntropy";
    E.channels  = cfg.channel;
    E.trials    = cfg.trials;
    E.times     = cfg.toi;
    E.freqs     = Freq;

    E.Energy  = Energy; % 每个时刻、每个尺度的能量
    E.Entropy = Entropy; 
    
    if(cfg.visualize)

    end

end

function [E] = getTotalWaveletEntropy(cfg,FT)

    N_trial   = numel(FT.trial);
    N_channel = numel(FT.label);
    N_order = cfg.order;

    % 离散小波变换
    % 注意，DWT是霍夫曼二叉树，从第n层的approx分出n+1层的approx和detail
    % 因此总能量=各层的detail+最后一层的approx
    Energy = cell(1,N_trial);
    Entropy = zeros(N_channel,N_trial);
    % 此外，需要注意：分解是从高频开始，低通/高通的分界线每次/=2
    % 第一次分解是Fs/4~Fs/2为高频、作为detail，0~Fs/4为低频
    ft_progress('init','etf');
    for i=1:N_trial
        N_time = size(FT.trial{i},2);

        ft_progress(i/N_trial,'正在处理试次(%d/%d)', i, N_trial);
        
        for j=1:N_channel
            Pow = zeros(N_time,N_order+1);
            Freq = zeros(1,N_order+2); % 因为每次/=2，所以这里存储频段的边界而非中心

            [C,L] = wavedec(FT.trial{i}(j,:),N_order,cfg.wavelet);
            
            % 分解时是从高频到低频，但我们存储时按照从低频到高频排列
            Freq(1) = FT.fsample / 2;
            for k=1:N_order
                idx = N_order+2-k; % 要存储的位置
                Pow(:,idx) = wrcoef('d',C,L,cfg.wavelet,k).^2; % 重建小波系数，然后计算功率

                Freq(k+1) = Freq(k) / 2;
            end
            Freq = flip(Freq); % 翻转，变成从低频到高频
            Pow(:,1) = wrcoef('a',C,L,cfg.wavelet,N_order).^2; % (times,bands)

            % 计算能量
            Pow_scale = sum(Pow,1); % 消除时间，保留尺度
            Pow_total = sum(Pow_scale); % 全部消除
            Energy{i}(:,:,j) = Pow; % (times,bands,channels)，bands=orders+1
            
            % Relative Wavelet Energy=特定尺度下所有时刻的能量之和/总能量
            P = Pow_scale ./ Pow_total;

            % Total Wavelet Entropy，基于不同尺度的相对能量计算
            Entropy(j,i) = -sum(P.*log(P));
        end

    end
    ft_progress('close');

    % 输出 ================================================================

    E = [];
    E.method    = "TotalWaveletEntropy";
    E.channels  = cfg.channel;
    E.trials    = cfg.trials;
    E.freqs     = Freq;

    E.Energy  = Energy; % 每个时刻、每个尺度的能量
    E.Entropy = Entropy; 
    
    if(cfg.visualize)
        figure("Name","TotalWaveletEntropy");
        imagesc(E.trials,E.channels,E.Entropy);
        colorbar;
        clim([0 max(E.Entropy,[],"all")]);
        colormap(hot);
    end

end




function [E] = getTimeEvolutionTotalWaveletEntropy(cfg,FT)

    N_trial   = numel(FT.trial);
    N_channel = numel(FT.label);
    N_order = cfg.order;
    N_window = numel(cfg.toi);

    % 离散小波变换
    % 注意，DWT是霍夫曼二叉树，从第n层的approx分出n+1层的approx和detail
    % 因此总能量=各层的detail+最后一层的approx
    Energy = cell(1,N_trial);
    Entropy = zeros(N_window,N_channel,N_trial);
    % 此外，需要注意：分解是从高频开始，低通/高通的分界线每次/=2
    % 第一次分解是Fs/4~Fs/2为高频、作为detail，0~Fs/4为低频
    ft_progress('init','etf');
    for i=1:N_trial
        N_time = size(FT.trial{i},2);

        ft_progress(i/N_trial,'正在处理试次(%d/%d)', i, N_trial);
            
        for j=1:N_channel
            Pow = zeros(N_time,N_order+1);
            Freq = zeros(1,N_order+2); % 因为每次/=2，所以这里存储频段的边界而非中心

            [C,L] = wavedec(FT.trial{i}(j,:),N_order,cfg.wavelet);
            
            % 分解时是从高频到低频，但我们存储时按照从低频到高频排列
            Freq(1) = FT.fsample / 2;
            for k=1:N_order
                idx = N_order+2-k; % 要存储的位置
                Pow(:,idx) = wrcoef('d',C,L,cfg.wavelet,k).^2; % 重建小波系数，然后计算功率

                Freq(k+1) = Freq(k) / 2;
            end
            Freq = flip(Freq); % 翻转，变成从低频到高频
            Pow(:,1) = wrcoef('a',C,L,cfg.wavelet,N_order).^2; % (times,bands)
            
            % 逐个计算时间窗
            for n=1:N_window
                Idx_win = FT.time{i}>=cfg.toi(n)-cfg.twin & FT.time{i}<=cfg.toi(n)+cfg.twin;
    
                % 计算能量
                Pow_scale = sum(Pow(Idx_win,:),1); % 消除时间，保留尺度
                Pow_total = sum(Pow_scale); % 全部消除
                Energy{i}(n,:,j) = Pow_scale; % (windows,bands,channels)，bands=orders+1
                
                % Relative Wavelet Energy=特定尺度下所有时刻的能量之和/总能量
                P = Pow_scale ./ Pow_total;
    
                % Total Wavelet Entropy，基于不同尺度的相对能量计算
                Entropy(n,j,i) = -sum(P.*log(P));
            end
        end
        
    end
    ft_progress('close');

    % 输出 ================================================================

    E = [];
    E.method    = "TimeEvolutionTotalWaveletEntropy";
    E.channels  = cfg.channel;
    E.trials    = cfg.trials;
    E.freqs     = Freq;
    E.times     = cfg.toi;

    E.Energy  = Energy; % 每个时刻、每个尺度的能量
    E.Entropy = Entropy; 
    
    if(cfg.visualize)

    end

end
