%{
  --------------------- �������� ------------------------------------------

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function       settings = parameter_settings()

% 
settings.c           = 3e8;

% ���β���
settings.fc          = 3e9;                                      % �ز�Ƶ��
settings.lambda      = settings.c/settings.fc;                   % ����
settings.PRF         = 1e3;                                      % PRF
settings.PRI         = 1/settings.PRF;                           % PRI
settings.B           = 20e6;                                     % ����
settings.tau         = 10e-6;                                    % ����
settings.mu          = settings.B/settings.tau;                  % ��Ƶ��
settings.BT          = settings.B*settings.tau;                  % ʱ������
settings.Pm          = 512;                                      % ������
settings.CPI         = settings.Pm*settings.PRI;                 % ��ش���ʱ��

% ϵͳ����
settings.fs          = 2*settings.B;                             % ����Ƶ��
settings.ts          = 1/settings.fs;                            % �������
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