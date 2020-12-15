%{
  ------------ ʹ��sinc��ֵ�ķ������һ��Keystone�任 ---------------------
  ���룺Sp_tf    --- Nums*Numf, NumsΪ��ʱ���������
  �����Src_tf   --- Nums*Numf, NumfΪ��ʱ���������
  ��Բ�ͬ��tm_model��Ҫ��ȡ��ͬ�Ĳ�ֵ��ʽ���ܵõ���ȷ�Ľ��

  -------------------------------------------------------------------------
  [1] κҫ������״�����˶�Ŀ������������о����Ͼ���ѧ��2013
%}
function     Src_tf = KT_sinc(Sp_tf)
% ȫ�ֱ���
global  settings;
% ��ʼ��
Numf    = settings.Numf;
Nums    = settings.Pm;
Src_tf  = zeros(Nums,Numf);
% ��ʱ��Ƶ��
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