%{
  --------------------- 参数设置 ------------------------------------------

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function       settings = parameter_settings()

% 
settings.c           = 3e8;

% 波形参数
settings.fc          = 3e9;                                      % 载波频率
settings.lambda      = settings.c/settings.fc;                   % 波长
settings.PRF         = 1e3;                                      % PRF
settings.PRI         = 1/settings.PRF;                           % PRI
settings.B           = 20e6;                                     % 带宽
settings.tau         = 10e-6;                                    % 脉宽
settings.mu          = settings.B/settings.tau;                  % 调频率
settings.BT          = settings.B*settings.tau;                  % 时宽带宽积
settings.Pm          = 512;                                      % 脉冲数
settings.CPI         = settings.Pm*settings.PRI;                 % 相关处理时间

% 系统参数
settings.fs          = 2*settings.B;                             % 采样频率
settings.ts          = 1/settings.fs;                            % 采样间隔
settings.samples     = ceil(settings.tau/settings.ts);
settings.Numf        = 2048;
settings.Nums        = settings.Pm;
settings.tm_Model    = 0;                                        % tm_Model
%--------------------------------------------------------------------------
if settings.tm_Model == 0
    settings.tm      = (0:settings.Nums-1);
else
    settings.tm      = -settings.Nums/2 : (settings.Nums/2 - 1);
end
settings.tm          = settings.tm.*settings.PRI;
%--------------------------------------------------------------------------
% SNR
settings.SNR         = -30;                                      % dB
settings.Noise       = 0;
settings.NoiseOnly   = 0;

return