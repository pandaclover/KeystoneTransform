%{
  ---------------- Ambiguity Search by Enrtopy ----------------------------
  % 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function    M = VambNumSearch_Entropy_SKT(Src_tf)

% 全局变量
global settings;

% 获得Numf
Numf     = settings.Numf;
% Nums     = settings.Pm;

% 快时间频率
fr       = (0:Numf-1).*(settings.fs/Numf) - settings.fs/2;

% 慢时间
tm       = settings.tm;

% 盲速度
Vamb     = settings.lambda*settings.PRF/2;

% 用于构造速度模糊函数匹配矩阵
[Fr, Tm] = meshgrid(fr,tm);

% 这里不需要对虚拟慢时间进行计算

%------------------------速度模糊校正--------------------------------------
% 1) 确定速度模糊数搜索的范围
M_l      = -20;
M_r      = 20;
M_scope  = M_l:M_r;
E_kamb   = zeros(1,length(M_scope));                         % 初始化图像熵

% 2）开始对速度模糊数进行搜索
for ii = 1:length(M_scope)
    
    m     = M_scope(ii);
    
    % 构建匹配滤波器进行搜索  --- Nums*Numf
    h_com = exp(1i*4*pi*m*Vamb.*Fr.*Tm/settings.c) ...
          .* exp(-1i*4*pi*m*Vamb.*Fr.^2.*Tm/settings.fc/settings.c);    
      
    % 3) 对快时间频域-慢时间的信号进行补偿，然后进行RCMC/integration
    xc_tf = Src_tf.*h_com;
    
    xc_tt = ifft(xc_tf,[],2);
    
    % 将脉压结果取模并累加起来 ---- 1*Numf
    xc_tt_sum  = sum(abs(xc_tt),1);
    
    % 4）计算xc_tt的熵
    Value1     = abs(xc_tt_sum).^2;
    
    Value2     = sum(Value1);
    
    E_kamb(ii) = -sum((Value1/Value2).*log(Value1/Value2));
    
end % for ii = 1:length(M_scope)

% 5）确定速度模糊数
E_kamb         = 1./E_kamb;
[~,label]      = find(E_kamb == max(E_kamb));
M              = M_scope(label);

%{

figure(101)
plot(M_scope, E_kamb);
grid on
xlabel('速度模糊数');
ylabel('熵幅度');

%}

return