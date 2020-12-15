%{
  ------------ 使用sinc插值的方法完成一阶Keystone变换 ---------------------
  输入：Sp_tf    --- Nums*Numf, Nums为慢时间采样点数
  输出：Src_tf   --- Nums*Numf, Numf为快时间采样点数
  针对不同的tm_model需要采取不同的插值公式才能得到正确的结果

  -------------------------------------------------------------------------
  [1] 魏耀，宽带雷达高速运动目标检测与成像处理研究，南京大学，2013
%}
function     Src_tf = KT_sinc(Sp_tf)
% 全局变量
global  settings;
% 初始化
Numf    = settings.Numf;
Nums    = settings.Pm;
Src_tf  = zeros(Nums,Numf);
% 快时间频率
fr      = (0:Numf-1).*(settings.fs/Numf) - settings.fs/2;

%------------------------TM Model------------------------------------------
if settings.tm_Model == 0
    
    %----------------- settings.tm_Model == 0 -----------------------------
    row_vec = 0:Nums-1;
    
    for n = 1:Numf
        
        alpha = settings.fc/(fr(n) + settings.fc); 
        
        for m = 1:Nums
            
            Src_tf(m,n) = sinc(alpha*(m-1) - row_vec)*Sp_tf(:,n);
        
        end % for m = 1:Nums
    
    end % for n = 1:Numf
    %----------------------------------------------------------------------
    
else
    
    %----------------- settings.tm_Model == 1 -----------------------------
    row_vec = -Nums/2:Nums/2-1;
    
    for n = 1:Numf
        
        alpha = settings.fc/(fr(n) + settings.fc); 
        
        for m = -Nums/2:Nums/2-1
            
            Src_tf(m+Nums/2+1,n) = ...
                      sinc(alpha*m - row_vec)*Sp_tf(:,n);
                  
        end % for m = 1:Nums
    
    end % for n = 1:Numf
    %----------------------------------------------------------------------
    
end % if settings.tm_Model == 0

return