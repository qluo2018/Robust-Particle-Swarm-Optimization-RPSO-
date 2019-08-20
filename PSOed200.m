%Author: Luo Qiang
%Abstract: minimizing function by PSO
%version: 1.0.1
%date: 7 June, 2007
%输入： 
%   目标函数描述： ObjPara = struct('NCBP',[min, max], 'Dim', dim_objfunction,'Name', objfunctionname)
%   可选函数： 'Schaffer','Rastrigian','Rossenbrock','Six-humpCamelback',
%   'Sphere','2n-minima','Griewank'.在myFun.m中定义。
%    算法参数描述： AlgPara = struct('MaxIte',max_iteration, 'WeightType', flag_weight, 'MatchType', flag_gamatch,
%                               'IniRange',[ini_low, ini_high], 'Algorithm', flag_algorithm, 'Size', size_particle);
%       已经实现的算法有： SPSO, LPSO, FPSO, RPSO
%       最大迭代次数max_iteration,
%       所有粒子采用单一权重:flag_weight = 1 或者每个粒子不同的权重flag_weight = 2; 
%       flag_weightcontrol  采用什么方法控制惯性权重: LPSO,FPSO,RPSO
%       对于RPSO还要考虑基因与粒子的匹配方法：如果需要进行随机匹配: flag_gamatch = 1; 
%           适度度匹配 = 0; 轮盘式选择匹配 = 2。   
%       初始化范围： ini_low_particle, ini_high_particle 
%
%
%
function Find = PSOed200(ObjPara, AlgPara)

%%%basic PSO parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim_objfunction = ObjPara.Dim;
NCBP_min = ObjPara.NCBP(1);
NCBP_max = ObjPara.NCBP(2);
name_objfunction = ObjPara.Name;

max_iteration = AlgPara.MaxIte;
flag_weight = AlgPara.WeightType;
flag_algorithm = AlgPara.Algorithm;
flag_gamatch = AlgPara.MatchType;
ini_low_particle =  AlgPara.IniRange(1);
ini_high_particle = AlgPara.IniRange(2);
size_particle = AlgPara.Size;



trueBest = 0;
errorBand = 1e-4; %保留4位有效数字


%每一维上的X允许范围
abs_max_particle = 10;
%每一维上的速度的允许范围

abs_max_V = 10;
%每一维上的速度的初始化范围
ini_low_V = 0;
ini_high_V = 10;
%权重的初始化种子
ini_inertia_weight = 0.9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%flag_weight = 'LPSO';  %'FPSO', 'FCAPSO'
%***********************
iteration = 1;
%%%Initialization of PSO parameters


clear c1_local;  clear Particles;
clear c2_global; clear V; clear W;
clear F; clear LBest; clear GBest;
clear GACourse;clear iniPop;
 
%Initialize iniPop between -0.04 and +0.04
iniPop = rand(20,5)*0.08-0.04;
OldPop = iniPop;
%Siginificance of local information and global information
c1_local(:,iteration) = 2 * ones(size_particle,1);    
c2_global(:,iteration) = 2 * ones(size_particle,1);    

%Initialization of positions of agents
% agents are initialized between ini_low_particle and 
% ini_high_particle randomly
Particles(:,:,iteration) = ini_low_particle ...
    + (ini_high_particle - ini_low_particle) ...
    * rand(size_particle, dim_objfunction, 1);

%Initialization of velocities of agents
%Between ini_low_V , ini_high_V, (which can also be started from zero)
V(:,:,:) = ini_low_V + (ini_high_V - ini_low_V) ...
    * rand(size_particle,dim_objfunction,1);

%Initialization of inertia weight of agents
if (flag_weight == 1)
    %for uniW %%%%%%%%%%%%%%%%%%%%%  UNIFORM WEIGHT
    W(:,:,:) = ini_inertia_weight ...
    * ones(size_particle,dim_objfunction,1);
else
    %for divW %%%%%%%%%%%%%%%%%%%%%  DIVERSITY WEIGHT
    %每个粒子，不同，但是每个粒子的每一维是相同的
    ini_inertia_weight_high = 0.90;
    ini_inertia_weight_low = 0.40;
    W(:,:,:) = rand(size_particle,1) * ones(1,dim_objfunction)...
        * (ini_inertia_weight_high ...
        - ini_inertia_weight_low) + ini_inertia_weight_low; 
    ini_W = W(:,:,1);
end

%Compute the ObjFunction Value for each Particle
F(:,iteration) = ComputeObjFunctionForParticle(Particles(:,:,iteration), name_objfunction);

%每个粒子的局部最小值
LBest(:,:) = F(:,iteration);
%记录局部最小位置
local_minimum_particle(:,:,:) =  Particles(:,:,iteration);

%所有粒子的全局最小值
[GBest(iteration),GIndex(iteration)] = min(LBest(:,iteration));  
%记录全局最小位置为了更新方便，特意形成一个每一页上每一行都相同的三维矩阵
global_minimum_particle(:,:,:) = ones(size_particle,1) * Particles(GIndex(iteration),:,iteration);

%%记录GA的演化过程
GACourse(iteration) = 0;
%%记录GA的适应度
GAFitness(:,iteration) = zeros(1,size_particle);
%%记录在单一权重算法中，每次选中的规则索引
CG_Index = 0;
%%记录被fire规则的适应度函数值
Index_Fitness = 0;
%-----------------------------------------%



%-----------------------------------------%






%%画图
% hpso = figure('Name','PSO');
% if (dim_objfunction == 2) 
%     subplot(2,3,1)
%     ezcontour('x^2 - 10 * cos(2*pi*x)+ y^2 - 10 * cos(2*pi*y) + 20',[-5,5],[-5,5])
%     title('Particles Flying')
%     hold on
%     h = plot(Particles(:,1,iteration),Particles(:,2,iteration),'o');
%     set(h, 'LineWidth',2)
%     axis([-5 5 -5 5])                
% end
%-------------------------------------

%迭代过程
no_decrease_count = 0;
term_condition = 1;
while(term_condition)    

    %update velocity
    V(:,:,iteration+1) = W(:,:,iteration) .* V(:,:,iteration)...
        - rand * (c1_local(:,iteration) * ones(1,dim_objfunction)) ...
        .* ( Particles(:,:,iteration) - local_minimum_particle(:,:,iteration) ) ...
        - rand * (c2_global(:,iteration) * ones(1,dim_objfunction))...
        .* ( Particles(:,:,iteration) - global_minimum_particle(:,:,iteration) );     
    %越界检查V
    V(:,:,iteration+1) = min(V(:,:,iteration+1),abs_max_V);  %让全部速度小于abs_max_V
    V(:,:,iteration+1) = max(V(:,:,iteration+1),-abs_max_V); %让全部速度大于-abs_max_V
    
    %particle fly
    Particles(:,:,iteration+1) = Particles(:,:,iteration) + V(:,:,iteration+1);
    %越界检查X
    Particles(:,:,iteration+1) = min(Particles(:,:,iteration+1),abs_max_particle);
    Particles(:,:,iteration+1) = max(Particles(:,:,iteration+1),-abs_max_particle);
    
%     %每次迭代都加入一定数量的随机粒子，以保持多样性
%     Particles(1:5,:,iteration+1) = ini_low_particle ...
%     + (ini_high_particle - ini_low_particle) ...
%     * rand(5, dim_objfunction, 1);
    
    
    %Compute the ObjFunction Value for each Particle
    F(:,iteration+1) = ComputeObjFunctionForParticle(Particles(:,:,iteration+1), name_objfunction);
  
    %update local best 
    for i = 1:size_particle
        if ( F(i,iteration+1) < LBest(i,iteration) )
            %如果最优值小了，更新            
            LBest(i,iteration+1) = F(i,iteration+1);
            local_minimum_particle(i,:,iteration+1) = Particles(i,:,iteration+1);            
        else
            %如果最优值没有变小，不更新            
            LBest(i,iteration+1) = LBest(i,iteration);
            local_minimum_particle(i,:,iteration+1) = local_minimum_particle(i,:,iteration);            
        end
    end
    
    %update global best 
    [GBest(iteration+1),GIndex(iteration+1)] = min(LBest(:,iteration+1));  
    if ( GBest(iteration+1) < GBest(iteration) )
        %如果变小，更新
        global_minimum_particle(:,:,iteration+1) = ones(size_particle,1) * Particles(GIndex(iteration+1),:,iteration+1);        
    else
        %否则，保持原来的全局最小
        global_minimum_particle(:,:,iteration+1) = global_minimum_particle(:,:,iteration);
        GBest(iteration+1) = GBest(iteration);
    end
    %记录减小量并计数
    decrease_iteration(iteration) = GBest(iteration+1) - GBest(iteration);
    if (decrease_iteration < 1e-10 )
        no_decrease_count = no_decrease_count + 1;       
    end    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%--------------change the inertia weight------------------------------
    switch (flag_algorithm)
        case 'SPSO'   
    %%%%%%%
     %%%%%%%
     plotkey = 0;
     %%%%%%%
     %%%%%%%
            %% 恒定权重
            %Method Ia：SimplePSO with uniW
            W(:,:,iteration+1) = W(:,:,iteration);
            %Method Ib: SimplePSO with divW
            %%---------------------------------------------------------------------           
        case 'LPSO'
             %%%%%%%
     %%%%%%%
     plotkey = 0;
     %%%%%%%
     %%%%%%%
            %%线性递减权重   
            if (flag_weight == 1)
                %Method IIa：LinearPSO
                W(:,:,iteration+1) = uniUpdateWeight(iteration+1,max_iteration, size_particle,dim_objfunction); %LPSO     
            else
                %Method IIb：LinearPSO
                W(:,:,iteration+1) = divUpdateWeight(iteration+1,max_iteration, ini_W);
            end           
             %%------------------------------------------------------------
        case 'FPSO'
             %%%%%%%
     %%%%%%%
     plotkey = 0;
     %%%%%%%
     %%%%%%%
            %%模糊调整权重           
            if (flag_weight == 1)
                %Method IIIa: FuzzyPSO
                CBPE = [NCBP_min, NCBP_max, GBest(iteration+1)];   %预期函数最小值和允许函数最大值，当前的最优值
                W_next = uniFAUpdateWeight(CBPE, W(1,1,iteration)); %整个权重矩阵应该就是一个数
                W(:,:,iteration+1) = W_next * ones(size_particle, dim_objfunction, 1);
            else
                %Method IIIb:
                CBPE = [NCBP_min, NCBP_max, GBest(iteration+1)]; 
                W(:,:,iteration+1) = divFAUpdateWeight(CBPE, W(:,:,iteration));%整个权重矩阵每一列都相同，即每个粒子的每一维共用一个权重
            end            
            %%-----------------------------------------------------------
            
        case 'RPSO'
            %%%%%%%
            %%%%%%%
            plotkey = 1;
            %%%%%%%
            %%%%%%%

            %%模糊文化算法调整权重
            if (flag_weight == 1)
%               %Method IVa: RPSO            
            else
                %Method IVb: RPSO  
                clear Fit_paras;
                Fit_paras = struct('sNCBP',[NCBP_min,NCBP_max],'sGBest_Curr', GBest(iteration),...
                    'sSizePSO',[size_particle,dim_objfunction],...
                    'sc1_local_current',c1_local(:,iteration),'sc2_global_current',c2_global(:,iteration), ...
                    'sV_Curr', V(:,:,iteration), 'sParticles_Curr', Particles(:,:,iteration), ...
                    'slocal_minimum_particle_Curr', local_minimum_particle(:,:,iteration) , ...
                    'sglobal_minimum_particle_Curr', global_minimum_particle(:,:,iteration), ... 
                    'sabs_max_V', abs_max_V, 'sabs_max_particle',abs_max_particle,...
                    'sLBest_Curr', LBest(:,iteration),'sW_Current',W(:,:,iteration),'sName',name_objfunction);
                GA_paras = struct('snvarsnvars',5,'sLB',-0.04*ones(5,1),'sUB',0.04*ones(5,1),...
                    'spopsize',size_particle, 'smaxGene',50,'iniPop',iniPop);
                clear delta_g;                 
                %计算适应度函数
                if iteration == 1
                    %计算初始CG的适应度
                    %让每条规则都fire一次，由此来计算其适应度
                    GAFitness(:,1) = InitilizeFitness(iniPop,Fit_paras);
                else
                    %计算fire规则后函数值改变了多少
                    GAFitness(:,iteration) = LBest(:,iteration) - F(:,iteration+1);                                      
                end
                %正规化到[-min,max],采用线性适应度定标      
                %belta = 0.05; %+ 100 * iteration / max_iteration;
%                 minGAFitness = min( GAFitness(:,iteration));
%                 maxGAFitness = max( GAFitness(:,iteration));
%                 difference = maxGAFitness - minGAFitness;
%                 if difference < 0.01
%                     difference = 1;
%                 end
%                 GAFitness(:,iteration) =  exp( (GAFitness(:,iteration) - ...
%                    minGAFitness ) / difference * 10 );
              %    GAFitness(:,iteration) =  0.05 * GAFitness(:,iteration) + 1000;
                   GAFitness(:,iteration) =  exp(GAFitness(:,iteration)); %对于Griewank这样的函数，本来值就小

               
                                
                %根据当前的适应度函数，产生新种群newPop             
                if (flag_gamatch == 1)                    
                    %随机找两个CG执行遗传操作,indexMatrix的第几行对应了第几个粒子匹配的两个规则索引                    
                    [Ytemp, indexMatrix] = sort(rand(size_particle,2));
                    %生成新的种群
                    NewPop = GeneticOperation(OldPop, indexMatrix);
                elseif (flag_gamatch == 0)
                    %固定找相邻两个CG执行遗传操作
                    indexMatrix = [1:size_particle;2:size_particle,1]';
                    %生成新的种群
                    NewPop = GeneticOperation(OldPop, indexMatrix);                    
                else               
                    %轮盘选择两个CG执行遗传操作
                    indexMatrix = RouletteSelect(GAFitness(:,iteration));                    
                    %生成新的种群
                    NewPop = GeneticOperation(OldPop, indexMatrix);                
                end

                CBPE = [NCBP_min, NCBP_max, GBest(iteration+1)];
                %采用NewPop更新惯性权重,已经完成了匹配
                W(:,:,iteration+1) = divFAUpdateWeight_RPSO(CBPE, W(:,:,iteration), NewPop); 
            end
            %%-----------------------------------------------------------
            %%更新OldPop
            OldPop = NewPop;

        otherwise
            'wrong flag'
    end
   
    %change the c1,c2   
    c1_local(:,iteration+1) =  c1_local(:,iteration);    
    c2_global(:,iteration+1) =  c2_global(:,iteration);   
    
    
    iteration = iteration + 1;
    term_condition = (iteration < max_iteration);
             %  && (no_decrease_count < 300 )); 
             %  && ( abs(GBest(iteration-1) - trueBest) > errorBand )  ); 
             
    %%%——画图————————————————————
    %%%%%
    %计算最优的方向与真实运动方向的夹角,最优值是0
%     for i_th = 6:20
%         a_theta = Particles(i_th,:,iteration-1);
%         b_theta = Particles(i_th,:,iteration) - Particles(i_th,:,iteration-1);
%         
%         
%         
%         if( norm(b_theta) < 10e-2 )
%             if ( norm(a_theta) < 10e-2 )
%                 theta(i_th,iteration) = 0;
%             else
%                 theta(i_th,iteration) = pi;
%             end               
%         else
%             tempp = sum(a_theta .* b_theta) / (norm(a_theta) * norm(b_theta));
%             if tempp >= -1 && tempp <= 1
%                 theta(i_th,iteration) = acos (tempp);
%             else
%                 i_th
%             end
%             
%         end
%     end    
%     
%     %画出粒子飞行的情况   
%     figure(hpso);
%     pause(0.0001);
%     if (iteration < max_iteration)
%         if (dim_objfunction == 2)
%             set(h,'Visible','off');
%         end
%     end    
%     if (dim_objfunction == 2)
%         subplot(2,3,1)
%         h = plot(Particles(:,1,iteration),Particles(:,2,iteration),'o');
%         set(h, 'LineWidth',2)
%         axis([-5 5 -5 5])           
%     end   
%     subplot(2,3,2)
%     bar(global_minimum_particle(1,:,iteration));
%     title('Best Solution')
%     subplot(2,3,3)
%     pr = minmax(LBest');
%     errorbar(1:iteration, mean(LBest),pr(:,1),pr(:,2));
%     colormap hsv;
%     %axis([1 iteration 0 600])
%     title('min, max and mean of the Best fval')
%     subplot(2,3,4)
%     plot(GBest);
%     title(['Global Best fval ',num2str(GBest(iteration))]);
%     subplot(2,3,5)
%     clear w_temp;
%     w_temp(:,:) = reshape(W(:,1,:),size_particle,iteration)';
%    %plot(1:iteration,w_temp(:,1),'g-',1:iteration, w_temp(:,2),'b-');
%     plot(1:iteration,w_temp(:,:));
%     title('inertia weight')
%     subplot(2,3,6)
%     plot(std(LBest));
%     title('std of the Best fval')
%     drawnow   
    % 录像使用
%     if (iteration == 2) 
%         pause(10); 
%     end
    %%%——————————————————————————   
        
end

%%%-----画图,对偏角进行统计分析,平均偏角的变化
% figure('Name','Angle')
% plot(1:size(theta,2), mean(theta));
%%%-------------------------------------------

Find = struct('best',GBest(iteration-1),...
    'BestGA',max(GAFitness), 'AveGA', mean(GAFitness),...
    'BestPSO', GBest(1:iteration-1), 'AvePSO', mean(LBest(:,1:iteration-1))); 

