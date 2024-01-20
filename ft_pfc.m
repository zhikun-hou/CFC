function [CFC] = ft_pfc(cfg,FT)


    % 检查参数 ================

    FT  = ft_checkdata(FT,'datatype','raw','feedback','yes');

    % 处理 ================================================================
    
    cfg  = ft_checkopt(cfg,'toi','ascendingdoublevector');
            
    cfg = ft_checkconfig(cfg, 'required',  {'bandphase','bandfreq','toi'});
    
    cfg  = ft_checkopt(cfg,'bandphase', 'ascendingdoublebivector'); % 带通滤波的[下限, 上限]
    cfg  = ft_checkopt(cfg,'bandfreq', 'ascendingdoublebivector');
    cfg  = ft_checkopt(cfg,'toi','ascendingdoublebivector');

    if(nargout==0)
        getPFC(cfg,FT);
    else
        CFC = getPFC(cfg,FT);
    end

end


function [PFC] = getPFC(cfg,FT)
    
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