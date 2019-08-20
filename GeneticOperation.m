%%遗传操作，交叉，变异
%%输入：父代种群及选择索引矩阵，该矩阵的每一行表示，子代父母的索引
%%输出：新的种群
function NewPop = GeneticOperation(OldPop, indexMatrix)

[PopSize, Length] = size(OldPop);
Bound = [-0.04,0.04];  %变异的范围不超过

%交叉，支持多点交叉
num_point = 2;
clear father;
clear mather;

for i = 1:PopSize
    father = OldPop(indexMatrix(i,1),:);
    mather = OldPop(indexMatrix(i,2),:);
    
    %随机产生交叉点坐标
    ind_point = 2 + floor( rand(num_point) * (Length-2) );
    %进行交叉
    offspring = father;
    if floor(num_point /2) == 0
        %交叉点个数为偶数
        for j = 1:2:num_point
            %对交叉点1，2之间的交叉，3，4之间的交叉，。。。
            offspring = [offspring(1:ind_point(j)-1),mather(ind_point(j):ind_point(j+1)),...
                offspring(ind_point(j+1)+1:Length)];
        end
    else
        %交叉点个数为奇数
        for j = 1:2:num_point-1
            %对交叉点1，2之间的交叉，3，4之间的交叉，。。。
            offspring = [offspring(1:ind_point(j)-1),mather(ind_point(j):ind_point(j+1)),...
                offspring(ind_point(j+1)+1:Length)];
        end
        %最后一个交叉点到最后都交换
        offspring = [offspring(1:ind_point(num_point)-1),mather(ind_point(num_point):Length)];
    end
    
    %单点一致变异，浮点数编码，需要知道上下界[-0.04,0.04]
    u = 1 + floor( rand(1) * Length );
    offspring(u) = 0.08*rand(1) - 0.04; 
    
    %
    NewPop(i,:) = offspring;
end
            
            
