%%调整惯性权重的多种方法
%%FPSO: 权重由模糊控制器来调节
% 输入：
%    CBPE: 预期函数最小值和允许函数最大值，当前的全局最优值
%    Weight_Current：当前使用的weight矩阵，其每一行都不同，但是每一列是相同的，
%                    即每个粒子的所有维都共用一个权重
% 输出：
%    W_Next: weight矩阵，其每一行都不同对应不同粒子的权重
% 
function W_Next = divFAUpdateWeight_FCAPSO( CBPE, Weight_Current, delta_g )

[size_particle, dim_obj] = size(Weight_Current);
%%  模糊规则的参数
%%  变量          Low     Medium   High
%%  NCBPE_i      (0,g1)   (g2,g3) (g4,1)
%%  a_i(t)       (0.2,g5) (g6,g7) (g8,1.1)
%%  DELTAa_i(t)   g9       g10     g11

%%默认参数
g = [0.055,0.055,0.35,0.35,0.5,0.5,0.75,0.75,0,0,0];

g = g + delta_g;

%%参数设置
para_NCBPE = [0,g(1); g(2),g(3); g(4),1];
para_Weight = [0.2, g(5); g(6),g(7); g(8),1.1];


%%  约束条件
%%  0 < g2 < g1 且 g4 < g３ < １
%%  0.2  < g6 < g5 且 g8 < g7 < 1.1
%%  -0.12 <= g9 < g10 < g11 <= 0.12
%% 约束条件检查
con1 = (0<g(2)) && (g(2)<g(1)) && (g(4)<g(3)) && (g(3)<1);
con2 = (0.2  < g(6) ) && (g(6) < g(5)) && (g(8) < g(7) ) && (g(7) < 1.1);
con3 = ( -0.12 <= g(9) ) && (g(9) < g(10)) && (g(10) < g(11)) && (g(11) <= 0.12 );

if (~con1 || ~con2 || ~con3)
    'conditions not meet for the parameters of the fuzzy rule base'
    exit(1);
end


%%  规则库  3*3 = 9条规则
% RuleBase = [1,1,2; 2,1,3; 3,1,3; 1,2,1; 2,2,2; 3,2,2; 1,3,1; 2,3,1; 3,3,1];
%RuleBase = [2,2,1;3,2,1;3,2,2];%2007/3
RuleBase = [2,1,1;3,2,1;3,2,2];%2007/2
        
%计算NCBPE
NCBPE = (CBPE(3) - CBPE(1)) / (CBPE(2) - CBPE(1));

%首先计算所有隶属度
mu_NCBPE =  Membership_Function(NCBPE, para_NCBPE); 
for i = 1 : size_particle
    mu_Weight(i,:) = Membership_Function(Weight_Current(i,1),para_Weight);
end
mu_DeltaW = [g(9), g(10), g(11)];

%模糊控制器
%针对每一个粒子都有一个权重
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

%更新权重，首先针对第一列更新
Weight_Current(:,1) = Weight_Current(:,1) + Delta_Weight(:,1);
%然后构造整个权重矩阵

W_Next = Weight_Current(:,1) * ones(1, dim_obj);