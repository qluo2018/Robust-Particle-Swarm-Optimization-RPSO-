%%��������Ȩ�صĶ��ַ���
%%FPSO: Ȩ����ģ��������������
% ���룺
%    CBPE: Ԥ�ں�����Сֵ�������������ֵ����ǰ��ȫ������ֵ
%    Weight_Current����ǰʹ�õ�ͳһ��һ��weightֵ
% �����
%    W_Next: ��һ��ʹ��ͳһ��һ��weightֵ
% 
function W_Next = uniFAUpdateWeight_FCAPSO( CBPE, Weight_Current, delta_g )

%%  ģ������Ĳ���
%%  ����          Low     Medium   High
%%  NCBPE_i      (0,g1)   (g2,g3) (g4,1)
%%  a_i(t)       (0.2,g5) (g6,g7) (g8,1.1)
%%  DELTAa_i(t)   g9       g10     g11

%%Ĭ�ϲ���
%g = [0.06,0.05,0.4,0.3,0.6,0.4,0.9,0.6,-0.02,0,0.05]; 
%  ��ʼ��Χ��x1������0.055Ϊ���ĵ��� 
% %%   x2������0.35Ϊ���ĵ���   
% %%   x3������0.5Ϊ���ĵ���
% %%   x4������0.75Ϊ���ĵ���
% %%   x5������0Ϊ���ĵ���
g = [0.055,0.055,0.35,0.35,0.5,0.5,0.75,0.75,0,0,0];
g = g + delta_g; 


%%��������
para_NCBPE = [0,g(1); g(2),g(3); g(4),1];
para_Weight = [0.2, g(5); g(6),g(7); g(8),1.1];


%%  Լ������
%%  0 < g2 < g1  g4 < g�� < ��
%%  0.2  < g6 < g5  g8 < g7 < 1.1
%%  -0.12 <= g9 < g10 < g11 <= 0.12
%% Լ���������
con1 = (0<g(2)) && (g(2)<g(1)) && (g(4)<g(3)) && (g(3)<1);
con2 = (0.2  < g(6) ) && (g(6) < g(5)) && (g(8) < g(7) ) && (g(7) < 1.1);
con3 = ( -0.12 <= g(9) ) && (g(9) < g(10)) && (g(10) < g(11)) && (g(11) <= 0.12 );

if (~con1 || ~con2 || ~con3)
    'conditions not meet for the parameters of the fuzzy rule base'
    exit(1);
end


%%  �����  3*3 = 9������
% RuleBase = [1,1,2; 2,1,3; 3,1,3; 1,2,1; 2,2,2; 3,2,2; 1,3,1; 2,3,1; 3,3,1];
%RuleBase = [2,2,1;3,2,1;3,2,2];%2007/3
RuleBase = [2,1,1;3,2,1;3,2,2];%2007/2
        
%����NCBPE
NCBPE = (CBPE(3) - CBPE(1)) / (CBPE(2) - CBPE(1));

%���ȼ�������������
mu_NCBPE =  Membership_Function(NCBPE, para_NCBPE); 
mu_Weight = Membership_Function(Weight_Current,para_Weight);
mu_DeltaW = [g(9), g(10), g(11)];

%ģ��������
Delta_Weight = 0; 
sum_temp = 0;
for k = 1:3
    for l = 1:3
        temp = 0;
        temp = mu_NCBPE(k) * mu_Weight(l);
        Delta_Weight = Delta_Weight + temp * mu_DeltaW(RuleBase(k,l));     
        sum_temp = sum_temp + temp;
    end
end
Delta_Weight =  Delta_Weight / sum_temp;


%����Ȩ�أ�������Ե�һ�и���
W_Next = Weight_Current + Delta_Weight;