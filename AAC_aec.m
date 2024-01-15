clear;
clc;



% =========================================================================

% 使用FieldTrip生成模拟信号

cfg = [];
cfg.method     = 'amplow_amphigh'; % phalow_freqhigh % phalow_amphigh % amplow_amphigh
cfg.fsample    = 1000;
cfg.trllen     = 2;
cfg.numtrl     = 100;
cfg.output     = 'mixed';

cfg.s1.freq    = 7;
cfg.s1.phase   = 'random'; %phase differs over trials
cfg.s1.ampl    = 1;

cfg.s2.freq    = 30;
cfg.s2.phase   = 'random'; %phase differs over trials
cfg.s2.ampl    = 1;

cfg.noise.ampl = 0.1;

FT = ft_freqsimulation(cfg);


% =========================================================================

% 关心的频率范围
F = [1:50];

% 计算每个试次、每个频段的包络序列
for i=1:numel(FT.trial)
    i
    for j=1:numel(F)
        f = F(j);
        % 目标频率-2到目标频率+2带通滤波，随后hilbert接abs提取包络
        Env(:,j,i) = abs(hilbert(bandpass(FT.trial{i},[f-2,f+2],1000)));
    end
end

% 滑动窗口，对于窗口内任意时刻交叉的两个频率，计算包络在时刻间的相关性
T_window = 2; % 因为模拟的是稳定的信号，所以单个时间窗即可
N_window = T_window*1000;
N_step   = 10;
n = 1;
for t=1:N_step:numel(FT.time{1})-N_window+1
    t
    for i=1:numel(F)
        for j=1:numel(F)
            for k=1:numel(FT.trial)
                R(i,j,n,k) = corr(Env(t:t+N_window-1,i,k),Env(t:t+N_window-1,j,k)); % 公式3
            end
        end
    end
    n = n+1;
end

% Fisher's Z变换，试次间平均，逆变换
AAC = tanh(mean(atanh(R),4)); % 公式4


AAC2 = sign(AAC) .* AAC.^2; % 公式5

imagesc(F,F,AAC2); axis xy; colorbar;

