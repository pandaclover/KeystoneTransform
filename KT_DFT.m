%{
  ------------- 使用DFT+IFFT方法完成一阶Keystone变换 ----------------------
  输入：Sp_tf   --- Nums*Numf，Nums为慢时间采样点数
  输出：Src_tf  --- Nums*Numf，Numf为快时间采样点数
  
  --- 但是当tm_model == 1时还是会出现距离偏移的情况
  --- 根据多普勒频移的大小偏移的格点数分别为0，±1，±2
  -_-!!!
  -------------------------------------------------------------------------
  [1] 魏耀，宽带雷达高速运动目标检测与成像处理研究，南京大学，2013
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   Src_tf  = KT_DFT(Sp_tf)
% 全局变量
global  settings;

Numf     = settings.Numf;
Nums     = settings.Pm;

% 初始化
S_kf     = zeros(Nums,Numf);

% 快时间频率
fr       = (0:Numf-1).*(settings.fs/Numf) - settings.fs/2;

% Step1: DFT过程
alpha    = (settings.fc + fr)./settings.fc;

% 公式(3.16)中求和号的上下界
m_vec    = 0:Nums-1; 
  
%--------------------------------------------------------------------------
for frIndex = 1:Numf
    
    % 公式(3.16)中所定义的(1+\xi*frIndex)
    Xi = alpha(frIndex);
    
    for fdIndex = 1:Nums
        
        %------------------------------------------------------------------
        S_kf(fdIndex,frIndex) = ...
           exp(-1i*2*pi*Xi*(fdIndex-Nums/2-1).*m_vec./Nums)*Sp_tf(:,frIndex);
        % 频域的范围必须按照-PRF/2到PRF/2来计算才能得到正确的结果
        
        %------------------------------------------------------------------      
        
    end % for fdIndex = 1:Nums
    
end % for frIndex

% Step 2: IFT
Src_tf = ifft(S_kf,[],1);

return