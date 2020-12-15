%{
  --------------------- test KT -------------------------------------------
  ����״�ز���ѹ�ź�Ƶ�ף�
  [1] A Novel Method for Parameter Estimation of Space Moving Target
  [2] ����״�����˶�Ŀ������������о����Ͼ���ѧ��2013
  
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc; echo off;

% ȫ�ֱ���
global     settings;

% ��������
settings    = parameter_settings();
error_flag  = 0;

Numf        = settings.Numf;
Nums        = settings.Pm;
Samples     = settings.samples;
tp          = settings.tm;

Rbin        = settings.c/(2*settings.fs);
df          = 1/settings.CPI;
dv          = df*settings.lambda/2;

% �ٶ�ģ����
Vamb        = settings.PRF*settings.lambda/2;
% ���̽�����
Rmax        = (Numf - Samples)*Rbin;

%--------------------------------------------------------------------------
% Ŀ���������
Target.R   = 500*Rbin;

%--------------------------------------------------------------------------
% Target.V   = 0;                              % M = 0, fin = 0

% Target.V   = 200*dv - Vamb;                  % M = -1, fin = -390.6250
% Target.V   = 200*dv + Vamb;                  % M = 1, fin = -390.6250
% Target.V   = -200*dv - Vamb;                 % M = -1, fin = 390.6250
% Target.V   = -200*dv + Vamb;                 % M = 1, fin = 390.6250
% Target.V   = 236*dv - Vamb;                 % M = -1, fin = -460.9375

Target.V   = 200*dv - 5*Vamb;                % M = -5, fin = -390.6250
% Target.V   = 200*dv + 5*Vamb;                % M = 5, fin = -390.6250
% Target.V   = -200*dv - 5*Vamb;               % M = -5, fin = 390.6250
% Target.V   = -200*dv + 5*Vamb;               % M = 5, fin = 390.6250

% Target.V   = 200*dv;                         % M = 0, fin = -390.6250
% Target.V   = -200*dv;                        % M = 0, fin = 390.6250
% Target.V   = -236*dv;                        % M = 0, fin = 460.9375
% Target.V   = 236*dv;                         % M = 0, fin = -460.9375  
% Target.V   = -50*dv - Vamb;                  % M = 0, fin = 97.6563

% Target.V   = randi([-100,100])*10;             % �ٶ�������Χ��-1000��1000

%--------------------------------------------------------------------------

Target.M   = round(Target.V/Vamb);
Target.Rt  = Target.R + Target.V.*tp;

if max(Target.Rt) > Rmax || min(Target.Rt) <= 0
    disp('-_-');
    return;
end

%--------------------------------------------------------------------------
% ��ʱ��
Tf         = (0:Numf-1).*settings.ts - settings.tau/2;

% �ο��ź�
Xt         = exp(1i*pi*settings.mu.*Tf.^2) ...
           .*(-settings.tau/2<=Tf & Tf<settings.tau/2);

% �ز��ź�
Tau        = 2.*Target.Rt./settings.c;
Tr         = repmat(Tf,Nums,1) - Tau';
Phasepoint = pi*settings.mu.*Tr.^2 - 2*pi*settings.fc.*Tau';
Xr         = exp(1i*Phasepoint).*(-settings.tau/2<=Tr & Tr<settings.tau/2);

% �źŹ���
SigPower   = sum(sum(abs(Xr).^2))/numel(Xr);

% �������� --- ����
sigma      = SigPower/(10^(settings.SNR/10));

%--------------------------------------------------------------------------
% ����ѹ��
St_tf      = fftshift(fft(Xt,[],2),2);
Sr_tf      = fftshift(fft(Xr,[],2),2);
Sp_tf      = Sr_tf.*conj(repmat(St_tf,Nums,1));
Sp_tt      = ifft(Sp_tf,[],2);

%--------------------------------------------------------------------------

%% RCMC
tic;
% 1) Sinc
% Src_tf     = KT_sinc(Sp_tf);
% ע��sinc��ֵ�Ľ����Ҫfftshift

% 2) DFT + IFT
% Src_tf     = KT_DFT(Sp_tf);

% 3) ChirpZ
Src_tf    = KT_CZT(Sp_tf);

toc;

%% Ambiguity Search
M = VambNumSearch_Entropy_SKT(Src_tf);

if abs(M - Target.M) >= 1
    
    disp('Error Model I');
    
end

fr         = (0:Numf-1).*(settings.fs/Numf) - settings.fs/2;
[Fr, Tp]   = meshgrid(fr,tp);

H_COM      = exp(1i*4*pi*M*Vamb.*Fr.*Tp/settings.c) ...
           .* exp(-1i*4*pi*M*Vamb.*Fr.^2.*Tp/settings.fc/settings.c);

Scv_tf     = Src_tf.*H_COM;       
       
Scv_tt     = ifft(Scv_tf,[],2);

% Scv_ft     = fftshift(fft(Scv_tt,[],1),1);             % sinc��ֵ����
Scv_ft     = fft(Scv_tt,[],1);                         % Ƶ�򷽷�

[Index.FI,Index.RI] = find(abs(Scv_ft) == max(max(abs(Scv_ft))));

fd         = (Index.FI - 1)*df - settings.PRF/2;
vi         = -settings.lambda*fd/2;
V          = vi + M*Vamb;
Index.V    = V;
R          = (Index.RI - 1)*Rbin;
Index.R    = R;
Index.M    = M;

% msgbox('Message Box','�������','warn');
disp('�������');




