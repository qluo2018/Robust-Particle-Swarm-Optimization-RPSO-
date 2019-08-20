%%��������Ȩ�صĶ��ַ���
%%FPSO: Ȩ����ģ��������������
% ���룺
%    CBPE: Ԥ�ں�����Сֵ�����������ֵ����ǰ��ȫ������ֵ
%    Weight_Current����ǰʹ�õ�weight������ÿһ�ж���ͬ������ÿһ������ͬ�ģ�
%                    ��ÿ�����ӵ�����ά������һ��Ȩ��
% �����
%    W_Next: weight������ÿһ�ж���ͬ��Ӧ��ͬ���ӵ�Ȩ��
% 
function W_Next = divFAUpdateWeight_FCAPSO( CBPE, Weight_Current, delta_g )

[size_particle, dim_obj] = size(Weight_Current);
%%  ģ������Ĳ���
%%  ����          Low     Medium   High
%%  NCBPE_i      (0,g1)   (g2,g3) (g4,1)
%%  a_i(t)       (0.2,g5) (g6,g7) (g8,1.1)
%%  DELTAa_i(t)   g9       g10     g11

%%Ĭ�ϲ���
g = [0.055,0.055,0.35,0.35,0.5,0.5,0.75,0.75,0,0,0];

g = g + delta_g;

%%��������
para_NCBPE = [0,g(1); g(2),g(3); g(4),1];
para_Weight = [0.2, g(5); g(6),g(7); g(8),1.1];


%%  Լ������
%%  0 < g2 < g1 �� g4 < g�� < ��
%%  0.2  < g6 < g5 �� g8 < g7 < 1.1
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
for i = 1 : size_particle
    mu_Weight(i,:) = Membership_Function(Weight_Current(i,1),para_Weight);
end
mu_DeltaW = [g(9), g(10), g(11)];

%ģ��������
%���ÿһ�����Ӷ���һ��Ȩ��
Delta_Weight = [0]; 
for i = 1 : size_particle
    Delta_Weight(i,1) = 0; 
    sum_temp = 0;
    for k = 1:3
        for l = 1:3
            temp = 0;
            temp = mu_NCBPE(k) * mu_Weight(i,l);
            Delta_Weight(i,1) = Delta_Weight(i,1) + temp * mu_DeltaW(RuleBase(k,l));     
            sum_temp = sum_temp + temp;
        end
    end
    Delta_Weight(i,1) =  Delta_Weight(i,1) / sum_temp;
end

%����Ȩ�أ�������Ե�һ�и���
Weight_Current(:,1) = Weight_Current(:,1) + Delta_Weight(:,1);
%Ȼ��������Ȩ�ؾ���

W_Next = Weight_Current(:,1) * ones(1, dim_obj);