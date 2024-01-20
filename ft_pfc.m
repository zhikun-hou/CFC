function [CFC] = ft_pfc(cfg,FT)


    % 检查参数 ================

    FT  = ft_checkdata(FT,'datatype','raw','feedback','yes');

    % 处理 ================================================================
    
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
            
            cfg = ft_checkconfig(cfg, 'required',  {'bandphase','bandfreq','toi'});
            
            cfg  = ft_checkopt(cfg,'bandphase', 'ascendingdoublebivector'); % 带通滤波的[下限, 上限]
            cfg  = ft_checkopt(cfg,'bandfreq', 'ascendingdoublebivector');
            cfg  = ft_checkopt(cfg,'toi','ascendingdoublebivector');

            if(nargout==0)
                getPFC_valid(cfg,FT);
            else
                CFC = getPFC_valid(cfg,FT);
            end
    end

end


function [PFC] = getPFC_valid(cfg,FT)
    
    cfg_get         = [];
    cfg_get.latency = cfg.toi;

    FT = ft_selectdata(cfg_get,FT);

    % 提取相位
    cfg_1 = [];
    cfg_1.channel  = 'all';
    cfg_1.bpfilter = 'yes';
    cfg_1.bpfreq   = cfg.bandfreq;
    cfg_1.hilbert  = 'unwrap_angle';
    Phi = ft_preprocessing(cfg_1,FT);
    
    % 提取瞬时相位变化率
    cfg_2 = [];
    cfg_2.channel = 'all';
    cfg_2.absdiff = 'yes';
    dPhi = ft_preprocessing(cfg_2,Phi);
    
    % 评估瞬时相位变化率的周期性
    cfg_3 = [];
    cfg_3.method     = 'mtmfft';
    cfg_3.output     = 'pow';
    cfg_3.taper      = 'hanning';
    cfg_3.foilim     = cfg.bandphase;
    cfg_3.channel    = 'all';
    cfg_3.keeptrials = 'no';
    
    PFC = ft_freqanalysis(cfg_3, dPhi);
    

    if(nargout==0)
        figure("Name","Phase-Freq Coupling");
        subplot(2,2,1);
        plot(PFC.freq,PFC.powspctrm');
        title("dPhi功率谱 - 线性");
        subplot(2,2,2);
        loglog(PFC.freq,PFC.powspctrm');
        title("dPhi功率谱 - 对数坐标");
        title("dPhi曲线");
    end

end