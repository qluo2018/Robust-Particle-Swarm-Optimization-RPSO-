%%�Ŵ����������棬����
%%���룺������Ⱥ��ѡ���������󣬸þ����ÿһ�б�ʾ���Ӵ���ĸ������
%%������µ���Ⱥ
function NewPop = GeneticOperation(OldPop, indexMatrix)

[PopSize, Length] = size(OldPop);
Bound = [-0.04,0.04];  %����ķ�Χ������

%���棬֧�ֶ�㽻��
num_point = 2;
clear father;
clear mather;

for i = 1:PopSize
    father = OldPop(indexMatrix(i,1),:);
    mather = OldPop(indexMatrix(i,2),:);
    
    %����������������
    ind_point = 2 + floor( rand(num_point) * (Length-2) );
    %���н���
    offspring = father;
    if floor(num_point /2) == 0
        %��������Ϊż��
        for j = 1:2:num_point
            %�Խ����1��2֮��Ľ��棬3��4֮��Ľ��棬������
            offspring = [offspring(1:ind_point(j)-1),mather(ind_point(j):ind_point(j+1)),...
                offspring(ind_point(j+1)+1:Length)];
        end
    else
        %��������Ϊ����
        for j = 1:2:num_point-1
            %�Խ����1��2֮��Ľ��棬3��4֮��Ľ��棬������
            offspring = [offspring(1:ind_point(j)-1),mather(ind_point(j):ind_point(j+1)),...
                offspring(ind_point(j+1)+1:Length)];
        end
        %���һ������㵽��󶼽���
        offspring = [offspring(1:ind_point(num_point)-1),mather(ind_point(num_point):Length)];
    end
    
    %����һ�±��죬���������룬��Ҫ֪�����½�[-0.04,0.04]
    u = 1 + floor( rand(1) * Length );
    offspring(u) = 0.08*rand(1) - 0.04; 
    
    %
    NewPop(i,:) = offspring;
end
            
            
