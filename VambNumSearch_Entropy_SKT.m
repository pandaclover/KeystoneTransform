%{
  ---------------- Ambiguity Search by Enrtopy ----------------------------
  % 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function    M = VambNumSearch_Entropy_SKT(Src_tf)

% ȫ�ֱ���
global settings;

% ���Numf
Numf     = settings.Numf;
% Nums     = settings.Pm;

% ��ʱ��Ƶ��
fr       = (0:Numf-1).*(settings.fs/Numf) - settings.fs/2;

% ��ʱ��
tm       = settings.tm;

% ä�ٶ�
Vamb     = settings.lambda*settings.PRF/2;

% ���ڹ����ٶ�ģ������ƥ�����
[Fr, Tm] = meshgrid(fr,tm);

% ���ﲻ��Ҫ��������ʱ����м���

%------------------------�ٶ�ģ��У��--------------------------------------
% 1) ȷ���ٶ�ģ���������ķ�Χ
M_l      = -20;
M_r      = 20;
M_scope  = M_l:M_r;
E_kamb   = zeros(1,length(M_scope));                         % ��ʼ��ͼ����

% 2����ʼ���ٶ�ģ������������
for ii = 1:length(M_scope)
    
    m     = M_scope(ii);
    
    % ����ƥ���˲�����������  --- Nums*Numf
    h_com = exp(1i*4*pi*m*Vamb.*Fr.*Tm/settings.c) ...
          .* exp(-1i*4*pi*m*Vamb.*Fr.^2.*Tm/settings.fc/settings.c);    
      
    % 3) �Կ�ʱ��Ƶ��-��ʱ����źŽ��в�����Ȼ�����RCMC/integration
    xc_tf = Src_tf.*h_com;
    
    xc_tt = ifft(xc_tf,[],2);
    
    % ����ѹ���ȡģ���ۼ����� ---- 1*Numf
    xc_tt_sum  = sum(abs(xc_tt),1);
    
    % 4������xc_tt����
    Value1     = abs(xc_tt_sum).^2;
    
    Value2     = sum(Value1);
    
    E_kamb(ii) = -sum((Value1/Value2).*log(Value1/Value2));
    
end % for ii = 1:length(M_scope)

% 5��ȷ���ٶ�ģ����
E_kamb         = 1./E_kamb;
[~,label]      = find(E_kamb == max(E_kamb));
M              = M_scope(label);

%{

figure(101)
plot(M_scope, E_kamb);
grid on
xlabel('�ٶ�ģ����');
ylabel('�ط���');

%}

return