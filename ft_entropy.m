function [Entropy] = ft_entropy(cfg,FT)

    % 预处理 ==============================================================

    FT  = ft_checkdata(FT,'datatype','raw','feedback','yes');
    
    cfg.trials  = ft_getopt(cfg,'trials',1:numel(FT.trial));
    cfg.channel = ft_getopt(cfg,'channel',1:numel(FT.label));
    cfg.latency = ft_getopt(cfg,'timerange','minperiod'); % 默认使用所有试次中最短的试次对所有试次进行截取
    
    FT = ft_selectdata(cfg,FT);

    % 检查公共设置 =========================================================

    cfg  = ft_checkopt(cfg,'method','char');

    cfg.embedding = ft_getopt(cfg,'embedding',3);
    cfg.scale = ft_getopt(cfg,'scale',1);
    cfg = ft_checkopt(cfg,'scale','ascendingdoublevector'); % 如果是向量，那就进行多尺度分析
    if(cfg.scale(1)<1)
        ft_error("多尺度分析时，尺度设置应当>=1");
    end
    if(any(mod(cfg.scale,1)~=0))
        ft_error("多尺度分析时，尺度应当为整数");
    end

    % 如果没有指定是否可视化，并且nargout为0，那就自动执行可视化
    cfg.visualize = ft_getopt(cfg,'visualize',nargout==0);

    % 执行 ================================================================
    
    N_scale = numel(cfg.scale);
    ft_progress('init','etf');
    for k=1:N_scale
        % 注意，scale是去掉中心后的窗口半长
        win = cfg.scale(k)*2+1;
        
        ft_progress(k/N_scale,'正在处理尺度(%d/%d)', k, N_scale);
        
        % 变换尺度 ------------------------------
        
        FT_temp = FT;
        if(win~=3)
            for i=1:numel(FT.trial)
                % movmean没有strade，只能用这个了
                FT_temp.trial{i} = matlab.tall.movingWindow(@mean,win,FT.trial{i}','Stride',win)';
                FT_temp.time{i}  = matlab.tall.movingWindow(@mean,win,FT.time{i}','Stride',win)';
            end
        end

        % 计算尺度下的熵 ------------------------

        switch(cfg.method)
            case {'approximate','Approximate'} % 近似熵
                E = getApproximateEntropy(cfg,FT_temp);
            case {'sample','Sample'} % 样本熵，对近似熵的改进
                E = getSampleEntropy(cfg,FT_temp);
            case {'fuzzy','Fuzzy'} % 模糊熵，对样本熵的改进
                E = getFuzzyEntropy(cfg,FT_temp);
            case {'permutation','Permutation'} % 排列熵
                E = getPermutationEntropy(cfg,FT_temp);
            otherwise
                ft_error("未知的熵计算方法");
        end

        E.scale = win;
        Entropy(k) = E;
    end

    ft_progress('close');

end

function [E] = getFuzzyEntropy(cfg,FT)

    N_trial        = numel(FT.trial);
    N_channel      = numel(FT.label);
    N_dimension    = cfg.embedding;

    E          = [];
    E.method   = "FuzzyEntropy";
    E.channels = 1:N_channel;
    E.trials   = 1:N_trial;

    % 论文推荐阈值为0.1-0.25
    cfg.fuzzy_n = ft_getopt(cfg,'fuzzy_n',2);
    cfg.fuzzy_r = ft_getopt(cfg,'fuzzy_r',0.2);

    % 计算 ==================================================

    Embedding = cell(N_trial,2);
    Entropy   = zeros(N_channel,N_trial);

    for i=1:N_trial
        Data = FT.trial{i}';
        
        % 将序列嵌入相空间 ----------------------

        N_embedding  = size(Data,1)-N_dimension;
        normalizer   = N_embedding * (N_embedding-1);

        Embedding0_trial = zeros(N_embedding-1,N_dimension,N_channel);
        Embedding1_trial = zeros(N_embedding-1,N_dimension+1,N_channel);
        for k=1:N_embedding-1
            Data0 = Data(k:k+N_dimension-1,:);
            Data1 = Data(k:k+N_dimension,:);
            Embedding0_trial(k,:,:) = Data0 - mean(Data0,1);
            Embedding1_trial(k,:,:) = Data1 - mean(Data1,1);
        end
        Embedding{i,1} = Embedding0_trial;
        Embedding{i,2} = Embedding1_trial;

        % 计算模糊熵 ----------------------------
            
        for j=1:N_channel
            X0 = Embedding0_trial(:,:,j);
            X1 = Embedding1_trial(:,:,j);
            
            % embedding两两之间计算切比雪夫距离
            % <=不可以放在squareform里，会有bug
            D0 = squareform(pdist(X0,'chebychev'));
            D1 = squareform(pdist(X1,'chebychev'));
            % 计算Fuzzy距离
            Mu0 = exp(-D0.^cfg.fuzzy_n ./cfg.fuzzy_r);
            Mu1 = exp(-D1.^cfg.fuzzy_n ./cfg.fuzzy_r);
            % 计算熵
            Phi0 = ( sum(Mu0,"all")-N_embedding ) / normalizer;
            Phi1 = ( sum(Mu1,"all")-N_embedding ) / normalizer;
            Entropy(j,i) = log(Phi0) - log(Phi1);
        end
    end

    % 输出 ===================================================

    E.Embedding = Embedding;
    E.Entropy   = Entropy;

    if(cfg.visualize)
        figure("Name","FuzzyEntropy");
        imagesc(E.trials,E.channels,E.Entropy);
        colormap(hot);
        colorbar;
    end

end


function [E] = getSampleEntropy(cfg,FT)

    N_trial        = numel(FT.trial);
    N_channel      = numel(FT.label);
    N_dimension    = cfg.embedding;

    E          = [];
    E.method   = "SampleEntropy";
    E.channels  = cfg.channel;
    E.trials    = cfg.trials;

    % 论文推荐阈值为0.1-0.25
    cfg.threshold = ft_getopt(cfg,'threshold',0.2);

    % 计算 ==================================================

    Embedding = cell(N_trial,1);
    Entropy   = zeros(N_channel,N_trial);

    for i=1:N_trial
        Data = FT.trial{i}';
        
        % 将序列嵌入相空间 ----------------------

        N_embedding     = size(Data,1)-N_dimension;
        Embedding_trial = zeros(N_embedding,N_dimension+1,N_channel);
        for k=1:N_embedding
            Embedding_trial(k,:,:) = Data(k:k+N_dimension,:);
        end
        Embedding{i} = Embedding_trial;

        % 计算样本熵 ----------------------------
        
        % 容限阈值
        Tolerance = cfg.threshold * std(Data,[],1); % (1,channels)
            
        for j=1:N_channel
            X0 = Embedding_trial(:,1:end-1,j);
            X1 = Embedding_trial(:,:,j);
            
            % embedding两两之间是否超过距离
            % <=不可以放在squareform里，会有bug
            D0 = squareform(pdist(X0,'chebychev'))<=Tolerance(j);
            D1 = squareform(pdist(X1,'chebychev'))<=Tolerance(j);
            % 超过阈值的计数
            C0 = sum(D0,"all")-N_embedding; % 矩阵运算时注意减去自己和自己
            C1 = sum(D1,"all")-N_embedding;
            % 计算熵
            Entropy(j,i) = -log(C1./C0);

        end
    end

    % 输出 ===================================================

    E.Embedding = Embedding;
    E.Entropy   = Entropy;

    if(cfg.visualize)
        figure("Name","SampleEntropy");
        imagesc(E.trials,E.channels,E.Entropy);
        colormap(hot);
        colorbar;
    end

end

function [E] = getApproximateEntropy(cfg,FT)

    N_trial        = numel(FT.trial);
    N_channel      = numel(FT.label);
    N_dimension    = cfg.embedding;

    E          = [];
    E.method   = "ApproximateEntropy";
    E.channels  = cfg.channel;
    E.trials    = cfg.trials;

    % 论文推荐阈值为0.1-0.25
    cfg.threshold = ft_getopt(cfg,'threshold',0.2);

    % 计算 ==================================================

    Embedding = cell(N_trial,2);
    Entropy   = zeros(N_channel,N_trial);

    for i=1:N_trial
        
        Data = FT.trial{i}';

        % 将序列嵌入相空间 ----------------------

        N_embedding     = size(Data,1)-N_dimension+1; % m对应的
        Embedding0_trial = zeros(N_embedding,N_dimension,N_channel);
        Embedding1_trial = zeros(N_embedding-1,N_dimension+1,N_channel);
        for k=1:N_embedding-1
            Embedding0_trial(k,:,:) = Data(k:k+N_dimension-1,:);
            Embedding1_trial(k,:,:) = Data(k:k+N_dimension,:);
        end
        Embedding0_trial(k+1,:,:) = Data(k+1:k+N_dimension,:);
        Embedding{i,1} = Embedding0_trial;
        Embedding{i,2} = Embedding1_trial;

        % 计算近似熵 ----------------------------
        
        Tolerance = cfg.threshold * std(Data,[],1);
            
        for j=1:N_channel
            X0 = Embedding0_trial(:,:,j);
            X1 = Embedding1_trial(:,:,j);
            
            % embedding两两之间是否超过距离
            % <=不可以放在squareform里，会有bug
            D0 = squareform(pdist(X0,'chebychev'))<=Tolerance(j);
            D1 = squareform(pdist(X1,'chebychev'))<=Tolerance(j);
            % 超过阈值的概率
            C0 = mean(D0);
            C1 = mean(D1);
            % 计算熵
            % 近似熵对于log=0的情况没有处理，所以结果为NaN很正常
            Phi0 = mean( log(nonzeros(C0)) );
            Phi1 = mean( log(nonzeros(C1)) );
            Entropy(j,i) = Phi0 - Phi1;
        end
    end

    % 输出 ===================================================

    E.Embedding = Embedding;
    E.Entropy   = Entropy;

    if(cfg.visualize)
        figure("Name","ApproximateEntropy");
        imagesc(E.trials,E.channels,E.Entropy);
        colormap(hot);
        colorbar;
    end

end

function [E] = getPermutationEntropy(cfg,FT)

    N_trial        = numel(FT.trial);
    N_channel      = numel(FT.label);
    N_dimension    = cfg.embedding;

    E           = [];
    E.method    = "PermutationEntropy";
    E.channels  = cfg.channel;
    E.trials    = cfg.trials;
    
    % 计算 ==================================================

    Embedding         = cell(N_trial,1);
    Entropy           = zeros(N_channel,N_trial);
    Entropy_normalize = zeros(N_channel,N_trial);

    for i=1:N_trial
        Data = FT.trial{i}';

        % 将序列嵌入相空间 ----------------------

        N_embedding     = size(Data,1)-N_dimension+1;
        Embedding_trial = zeros(N_embedding,N_dimension,N_channel);
        for k=1:N_embedding
            Embedding_trial(k,:,:) = Data(k:k+N_dimension-1,:);
        end
        Embedding{i} = Embedding_trial;

        % 计算排列熵 ----------------------------

        for j=1:N_channel
            Embedding_channel = Embedding_trial(:,:,j);
            
            % 计算排序模式
            [~,Sort_idx] = sort(Embedding_channel,2);
            % 统计排序模式
            [~,~,Pattern_idx] = unique(Sort_idx,'rows');
            Pattern_count = groupcounts(Pattern_idx);
            % 计算每种排序模式的出现概率
            P = Pattern_count./sum(Pattern_count);
            % 计算熵
            H0 = -sum(P.*log(P));
            H1 = H0 / log(numel(P));
            Entropy(j,i)           = H0;
            Entropy_normalize(j,i) = H1;
        end

    end

    % 输出 ===================================================

    E.Embedding         = Embedding;
    E.Entropy           = Entropy;
    E.Entropy_normalize = Entropy_normalize;

    if(cfg.visualize)
        figure("Name","PermutationEntropy");
        imagesc(E.trials,E.channels,E.Entropy_normalize);
        colormap(hot);
        colorbar;
    end
end

