
### 文件说明
- 非ft开头：原型算法
- ft开头：直接可与fieldtrip联合使用的函数

- ft_entropy.m
  在time domain计算entropy，包括Approximate/Sample/Fuzzy/Permutation Entropy
- ft_entropyTF.m
  基于ft_freqanalysis的计算结果（傅里叶或STFT）在频域/时频域计算Spectral Entropy
- ft_entropyWT.m
  基于离散小波变换计算Wavelet Entropy以及相关的衍生指标
  
- AAC_aec.m
  使用Amplitude Envelope Correlation方法估计Amplitude-Amplitude Coupling
  【参考文献】Amplitude envelope correlation detects coupling among incoherent brain signals
