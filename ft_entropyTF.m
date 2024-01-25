
% 在时频域计算熵，输入ft_freqanalysis的计算结果

function [Entropy] = ft_entropyTF(cfg,TF)

    % 预处理 ==============================================================

    TF  = ft_checkdata(TF,'datatype','freq','feedback','yes');
    
    % 检查公共设置 =========================================================

    % 如果没有指定是否可视化，并且nargout为0，那就自动执行可视化
    cfg.visualize = ft_getopt(cfg,'visualize',nargout==0);

    % 执行 ================================================================
    
    switch(TF.cfg.output)
        case {'pow','powandcsd'}
            
        case {'fourier'}
            TF.powspctrm = abs(TF.fourierspctrm).^2;
        otherwise
            error("功率谱熵的输入应当为ft_freqanalysis输出的output=pow/powandcsd/fourier");
    end

    Entropy = getSpectralEntropy(cfg,TF);

end



function [E] = getSpectralEntropy(cfg,TF)

    E = [];
    E.method   = "SpectralEntropy";
    E.trials   = 1:size(TF.powspctrm,1);
    E.channels = 1:size(TF.powspctrm,2);
    E.freqs    = TF.freq;
    E.times    = TF.time;
    
    % 计算 ================================================================
    
    % matlab有自带的pentropy(基于pspectrum)和wentropy，但是算出来的结果完全不同

    switch(TF.cfg.method)
        case 'mtmfft' % FFT的功率谱熵
            
            Pow = TF.powspctrm; % keeptrials则(trials,channels,freqs)，否则(channels,freqs)

            if(strcmp(TF.dimord,'rpt_chan_freq'))
                P = Pow ./ sum(Pow,3);
                H0 = -sum(P.*log2(P),3)'; % (channels,trials)
            elseif(strcmp(TF.dimord,'chan_freq'))
                P = Pow ./ sum(Pow,2);
                H0 = -sum(P.*log2(P),2); % (channels,1)
            end

            % matlab自带的pentropy使用的是log2，为了保持一致性，这里也使用log2
            H1 = H0 ./ log2(numel(TF.freq));

        case 'mtmconvol' % STFT的功率谱熵
                
            Pow = TF.powspctrm; % keeptrials则(trials,channels,freqs,times)，否则(channels,freqs,times)

            if(strcmp(TF.dimord,'rpt_chan_freq_time'))
                P = Pow ./ sum(Pow,3);
                H0 = -sum(P.*log2(P),3); % (trials,channels,1,times)
                H0 = permute(H0,[4,2,1,3]); % (times,channels,trials)
            elseif(strcmp(TF.dimord,'chan_freq_time'))
                P = Pow ./ sum(Pow,2);
                H0 = -sum(P.*log2(P),2); % (channels,1,times)
                H0 = permute(H0,[3,1,2]); % (times,channels)
            end

            H1 = H0 ./ log2(numel(TF.freq));

        otherwise
            ft_error("不支持的方法");
    end
    
    % 输出 ================================================================

    E.Entropy           = H0;
    E.Entropy_normalize = H1;

    if(cfg.visualize)


    end

end