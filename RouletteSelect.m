%%轮盘选择
%%输入：种群的适应度函数
%%输出：轮盘选择的父母个体索引矩阵
function indexMatrix = RouletteSelect(GAFitness)
size_particle = size(GAFitness,1);
Roulette = GAFitness / sum(GAFitness);
u = rand(size_particle,2);
for i_index = 1:size_particle
    i_sum(i_index) = sum(Roulette(1:i_index));
end
i_sum(size_particle) = 1;
for i_matrix = 1 : size_particle
    for j_matrix = 1:2
        if u(i_matrix,j_matrix) <= i_sum(1)
            indexMatrix(i_matrix,j_matrix) = 1;   
        else
            for i_index = 2:size_particle
                if u(i_matrix,j_matrix) > i_sum(i_index-1)...
                        && u(i_matrix,j_matrix) <= i_sum(i_index)
                    indexMatrix(i_matrix,j_matrix) = i_index;
                    break;
                end
            end
        end
    end
end