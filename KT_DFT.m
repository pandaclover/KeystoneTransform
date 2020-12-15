%{
  ------------- ʹ��DFT+IFFT�������һ��Keystone�任 ----------------------
  ���룺Sp_tf   --- Nums*Numf��NumsΪ��ʱ���������
  �����Src_tf  --- Nums*Numf��NumfΪ��ʱ���������
  
  --- ���ǵ�tm_model == 1ʱ���ǻ���־���ƫ�Ƶ����
  --- ���ݶ�����Ƶ�ƵĴ�Сƫ�Ƶĸ�����ֱ�Ϊ0����1����2
  -_-!!!
  -------------------------------------------------------------------------
  [1] κҫ������״�����˶�Ŀ������������о����Ͼ���ѧ��2013
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   Src_tf  = KT_DFT(Sp_tf)
% ȫ�ֱ���
global  settings;

Numf     = settings.Numf;
Nums     = settings.Pm;

% ��ʼ��
S_kf     = zeros(Nums,Numf);

% ��ʱ��Ƶ��
fr       = (0:Numf-1).*(settings.fs/Numf) - settings.fs/2;

% Step1: DFT����
alpha    = (settings.fc + fr)./settings.fc;

% ��ʽ(3.16)����ͺŵ����½�
m_vec    = 0:Nums-1; 
  
%--------------------------------------------------------------------------
for frIndex = 1:Numf
    
    % ��ʽ(3.16)���������(1+\xi*frIndex)
    Xi = alpha(frIndex);
    
    for fdIndex = 1:Nums
        
        %------------------------------------------------------------------
        S_kf(fdIndex,frIndex) = ...
           exp(-1i*2*pi*Xi*(fdIndex-Nums/2-1).*m_vec./Nums)*Sp_tf(:,frIndex);
        % Ƶ��ķ�Χ���밴��-PRF/2��PRF/2��������ܵõ���ȷ�Ľ��
        
        %------------------------------------------------------------------      
        
    end % for fdIndex = 1:Nums
    
end % for frIndex

% Step 2: IFT
Src_tf = ifft(S_kf,[],1);

return