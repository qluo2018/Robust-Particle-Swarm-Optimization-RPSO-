%适应度函数
%浮点数编码
%输入： 
% GA 是最小化

%Fit_paras =
%struct('sNCBP',[NCBP_min,NCBP_max],'sGBest_Curr',GBest(iteration),'sSizePSO',[size_particles,dim_objfunction],...
%      'sc1_local_current',c1_local_current(:,iteration),'sc2_global_current',c2_global_current(:,iteration), ...
%      'sV_Curr', V_Current(:,:,iteration), 'sParticles_Curr', Particles(:,:,iteration), ...
%      'slocal_minimum_particle_Curr', local_minimum_particle(:,:,iteration) , ...
%      'sglobal_minimum_particle_Curr', global_minimum_particle(:,:,iteration), ... 
%      'sabs_max_V', abs_max_V,'sabs_max_particle',abs_max_particle,...
%      'sLBest_Curr', LBest(:,iteration),'sW_Current',W(:,:,iteration))
%输出：result_Fitness，文化基因的适应度
%%%% %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%% 
%%%% 每个基因的适应度是对整个粒子群体的效果，最后选出最好的基因来做为下一代使用的规则参数
%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%% 

function result_Fitness = divmyFitness ( x, Fit_paras )

result_Fitness = 0;

delta_g = decode_myga(x);

num_var = size(delta_g,2);
if (num_var ~= 11)
    'error in CG'
    exit(1)
end


size_particle = Fit_paras.sSizePSO(1);
dim_objfunction = Fit_paras.sSizePSO(2);
%update W
 CBPE = [Fit_paras.sNCBP(1), Fit_paras.sNCBP(2), Fit_paras.sGBest_Curr]; %计算NCBPE时我们还是针对整个种群来计算，减少计算量，只是每个粒子的惯性权重不一样
 W_NextIter = divFAUpdateWeight_FCAPSO(CBPE, Fit_paras.sW_Current, delta_g); %整个权重矩阵每一列都相同
 
%update velocity
    V_Next(:,:) = W_NextIter .* Fit_paras.sV_Curr...
        - rand * (Fit_paras.sc1_local_current * ones(1,dim_objfunction)) ...
        .* ( Fit_paras.sParticles_Curr - Fit_paras.slocal_minimum_particle_Curr ) ...
        - rand * (Fit_paras.sc2_global_current * ones(1,dim_objfunction))...
        .* ( Fit_paras.sParticles_Curr - Fit_paras.sglobal_minimum_particle_Curr );     
    %越界检查V
    V_Next = min(V_Next,Fit_paras.sabs_max_V);  %让全部速度小于abs_max_V
    V_Next = max(V_Next,-Fit_paras.sabs_max_V); %让全部速度大于-abs_max_V
    
    %particle fly
    Particles_Next(:,:) =  Fit_paras.sParticles_Curr + V_Next;
    %越界检查X
    Particles_Next = min(Particles_Next, Fit_paras.sabs_max_particle);
    Particles_Next = max(Particles_Next, -Fit_paras.sabs_max_particle);
    
    
    
    %Compute the ObjFunction Value for each Particle
    F_Next = ComputeObjFunctionForParticle(Particles_Next);
    
    %用局部最优值更新值与局部最优比较，更新适应度
    result_Fitness = 0;
    for i = 1:size_particle
        if ( F_Next(i) < Fit_paras.sLBest_Curr(i) )
            %如果最优值小了，更新            
            result_Fitness = result_Fitness + 1 + W_NextIter(i,1);   
             %如果最优值小了，更新            
            LBest_Next(i) = F_Next(i);             
        else
            %如果最优值没有变小，不更新            
            %result_Fitness = result_Fitness + 0;    
             %如果最优值没有变小，不更新            
            LBest_Next(i) = Fit_paras.sLBest_Curr(i);            
        end
    end
    result_Fitness = result_Fitness / size_particle;   %文化基因对于种群适应度修改的平均效果
    
    %如果全局最优值获得更新，那么这个基因的适应度比较高 
    [GBest_Next,GIndex_Next] = min(abs(LBest_Next));  
    if ( GBest_Next < Fit_paras.sGBest_Curr )
        %如果变小，更新
        result_Fitness = 3;       
    end
   result_Fitness = 10 - result_Fitness;  
 
 

