% Amplitude-Envelope Coupling

function [CFC] = ft_aec(cfg,varargin)
    
    % 检查参数 ============================================================
    
    if(nargin>3)
        error("Too many data");
    elseif(nargin==2)
        FT_1  = varargin{1};
        FT_2  = FT_1;

        FT_1  = ft_checkdata(FT_1,'datatype','raw','feedback','yes');
    else
        FT_1  = varargin{1};
        FT_2  = varargin{2};

        FT_1  = ft_checkdata(FT_1,'datatype','raw','feedback','yes');
        FT_2  = ft_checkdata(FT_2,'datatype','raw','feedback','yes');

        % 检验两组数据是否能够一一对应
        if(FT_1.fsample ~= FT_2.fsample)
            error("FT_1与FT_2的fsample不同");
        end
        if(numel(FT_1.trial) ~= numel(FT_2.trial))
            error("FT_1与FT_2的trial数量不同");
        end
        if(numel(FT_1.label) ~= numel(FT_2.label))
            error("FT_1与FT_2的channel数量不同");
        end
        for i=1:numel(FT_1.trial)
            if(numel(FT_1.time{i}) ~= numel(FT_2.time{i}))
                error("试次"+i+"中，FT_1与FT_2的时间长度不同");
            end
        end
    end
    
    % 处理设置 =============================================================
    
    cfg  = ft_checkopt(cfg,'toi','ascendingdoublevector');

    % 模式：探索性分析explore/验证性分析valid
    cfg.mode = ft_getopt(cfg,'mode',"explore"); % ft_getopt的key必须是char，不能是str
    
    switch(cfg.mode)
        case {"explore","Explore"}
            disp("【AEC】探索性分析");
            disp("对于对应通道的所有时刻，计算各频段两两之间的试次间包络相关");
            
            cfg = ft_checkconfig(cfg, 'required',  {'foi','fwidth','toi'});

            cfg  = ft_checkopt(cfg,'toi','ascendingdoublebivector'); % 探索模式下，toi只能是一个range，否则矩阵维度太高了
            cfg  = ft_checkopt(cfg,'foi','ascendingdoublevector');
            cfg  = ft_checkopt(cfg,'fwidth','numericscalar');
            
            CFC = getAEC_explore(cfg,FT_1,FT_2);
            
        
        case {"valid","Valid"}
            disp("【AEC】验证性分析");
            disp("对于数据集的对应通道，计算逐时刻的试次间包络相关");
            

            cfg = ft_checkconfig(cfg, 'required',  {'bandlow','bandhigh','toi'});
            
            cfg  = ft_checkopt(cfg,'bandlow', 'ascendingdoublebivector'); % 带通滤波的[下限, 上限]
            cfg  = ft_checkopt(cfg,'bandhigh','ascendingdoublebivector');
            cfg  = ft_checkopt(cfg,'toi','ascendingdoublevector'); % 验证模式下，toi可以是向量

            CFC = getAEC_valid(cfg,FT_1,FT_2);
    end
    

end

function [CFC] = getAEC_explore(cfg,FT_1,FT_2)

    % 根据timerange截取数据 --------------------------------------------
    cfg_get         = [];
    cfg_get.latency = [min(cfg.toi) max(cfg.toi)];

    FT_1 = ft_selectdata(cfg_get,FT_1);
    FT_2 = ft_selectdata(cfg_get,FT_2);

    % -----------------------------------------------------------------

    % 注意：由于FT_1和FT_2可能是不同的数据，探索分析时矩阵可能是不对称的

    CFC = [];
    N_freq    = numel(cfg.foi);
    N_channel = numel(FT_1.label);
    N_trial   = numel(FT_1.trial);
    AEC       = zeros(N_freq,N_freq,N_channel);

    for i=1:N_freq
        Band_1 = cfg.foi(i) + [ -cfg.fwidth , cfg.fwidth ];

        for j=1:N_freq
            Band_2 = cfg.foi(j) + [ -cfg.fwidth , cfg.fwidth ];
            [Env_1,Env_2] = getEnvelope(FT_1,Band_1,FT_2,Band_2);

            % 计算相关 -----------------------------------------------------
            % 注意，探索模式下只有一个时间段
            R = zeros(N_channel,N_trial);
  
            for k=1:N_trial
                for w=1:N_channel % each channel
                    R(w,k) = corr(Env_1(w,:,k)',Env_2(w,:,k)');
                end
            end
        
            % Fisher's Z Transform ----------------------------------------
            
            % 变换-试次间平均-逆变换
            AEC_temp = tanh(mean(atanh(R),2));
            AEC_temp = AEC_temp.^2 .* sign(AEC_temp);
            AEC(i,j,:) = AEC_temp;
        end
    end
       
    % 包装 ================================================================

    CFC = [];
    CFC.freq    = cfg.foi;
    CFC.channel = 1:N_channel;
    CFC.label   = FT_1.label;
    CFC.spctrm  = AEC;

end

function [CFC] = getAEC_valid(cfg,FT_low,FT_high)

    % 根据timerange截取数据 ---------------------------------------
    cfg_get         = [];
    cfg_get.latency = [min(cfg.toi) max(cfg.toi)];

    FT_low  = ft_selectdata(cfg_get,FT_low);
    FT_high = ft_selectdata(cfg_get,FT_high);
    
    % 带通滤波+希尔伯特 --------------------------------------------

    [Env_low,Env_high] = getEnvelope(FT_low,cfg.bandlow,FT_high,cfg.bandhigh);

    % 计算包络幅度的相关性 ------------------------------------------

    % 计算相关
    N_trial   = numel(FT_low.trial);
    N_channel = numel(FT_low.label);
    R = zeros(N_channel,numel(cfg.toi)-1,N_trial);
    Times = FT_low.time{1};

    for j=1:numel(cfg.toi)-1
        T_idx = Times>cfg.toi(j) & Times<cfg.toi(j+1);

        for i=1:N_channel % each channel
            for k=1:N_trial
                R(i,j,k) = corr(Env_low(i,T_idx,k)',Env_high(i,T_idx,k)');
            end
        end
    end

    % Fisher's Z Transform ----------------------------------------
    
    % 变换-试次间平均-逆变换
    AEC = tanh(mean(atanh(R),3));
    AEC = AEC.^2 .* sign(AEC);

    % 包装结果 -----------------------------------------------------

    CFC = [];
    CFC.time    = (cfg.toi(1:end-1)+cfg.toi(2:end))/2;
    CFC.channel = 1:N_channel;
    CFC.label   = FT_low.label;
    CFC.spctrm  = AEC;

end

% 带通滤波+希尔伯特，之后返回包络幅度矩阵
function [Env_1,Env_2] = getEnvelope(FT_1,Band_1,FT_2,Band_2)

    cfg_bp = [];
    cfg_bp.bpfilter   = 'yes';
    cfg_bp.keeptrials = 'yes';
    cfg_bp.hilbert    = 'abs'; % 取包络

    cfg_bp.bpfreq = Band_1;
    FT_1 = ft_preprocessing(cfg_bp,FT_1);
    cfg_bp.bpfreq = Band_2;
    FT_2 = ft_preprocessing(cfg_bp,FT_2);

    % ======================================

    N_trial   = numel(FT_1.trial);
    N_channel = numel(FT_1.label);
    N_time    = numel(FT_1.time{1});

    Env_1 = zeros(N_channel,N_time,N_trial);
    Env_2 = zeros(N_channel,N_time,N_trial);

    for i=1:N_trial
        Env_1(:,:,i) = FT_1.trial{i};
        Env_2(:,:,i) = FT_2.trial{i};
    end

end