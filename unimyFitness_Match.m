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
%      'sLBest_Curr', LBest(:,iteration),'sW_Current',W(1,1,iteration))
%输出：result_Fitness，文化基因的适应度
%%%% %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%% 
%%%% 每个基因与粒子1对1随机匹配
%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%% 

function result_Fitness = unimyFitness_Match ( x, Fit_paras )

%%%%%%%%%%%%%%%% 两个全局变量用于匹配基因和粒子
global fitnesstimes;
global MatchGenesWithParticles;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result_Fitness = 0;

delta_g = decode_myga(x);

num_var = size(delta_g,2);
if (num_var ~= 11)
    'error in CG'
    exit(1)
end
size_particle = Fit_paras.sSizePSO(1);
dim_objfunction = Fit_paras.sSizePSO(2);


%%%%%%%%%%%%%%%%使用全局变量%%%%%%%%%%%%%%%%%%%%%%%
index_temp = MatchGenesWithParticles(fitnesstimes); 


%update W
 CBPE = [Fit_paras.sNCBP(1), Fit_paras.sNCBP(2), Fit_paras.sGBest_Curr];
 W_NextIter = uniFAUpdateWeight_FCAPSO(CBPE, Fit_paras.sW_Current, delta_g); %整个权重矩阵应该就是一个数
 
 clear V_Next;
 clear Particles_Next;
 %update velocity
 V_Next = W_NextIter * Fit_paras.sV_Curr(index_temp,:)...
        - rand * Fit_paras.sc1_local_current(1) ...
        * ( Fit_paras.sParticles_Curr(index_temp,:) - Fit_paras.slocal_minimum_particle_Curr(index_temp,:) ) ...
        - rand * Fit_paras.sc2_global_current(1)...
        * ( Fit_paras.sParticles_Curr(index_temp,:) - Fit_paras.sglobal_minimum_particle_Curr(index_temp,:) );     

    %越界检查V
    V_Next = min(V_Next,Fit_paras.sabs_max_V);  %让全部速度小于abs_max_V
    V_Next = max(V_Next,-Fit_paras.sabs_max_V); %让全部速度大于-abs_max_V
    
    %particle fly
    Particles_Next =  Fit_paras.sParticles_Curr (index_temp,:) + V_Next;
    %越界检查X
    Particles_Next = min(Particles_Next, Fit_paras.sabs_max_particle);
    Particles_Next = max(Particles_Next, -Fit_paras.sabs_max_particle);
    

    
    %Compute the ObjFunction Value for each Particle
    F_Next = ComputeObjFunctionForParticle(Particles_Next);
    

   
       
   result_Fitness = 0;
    %用局部最优值更新值与局部最优比较，更新适应度
        if ( F_Next < Fit_paras.sLBest_Curr(index_temp) )
            %如果最优值小了，更新            
            result_Fitness = 1 +  W_NextIter;                              
        end
    %与整体最优比较，更新适应度      
    if ( F_Next < Fit_paras.sGBest_Curr )
        %如果变小，更新
        result_Fitness = 3;       
    end
    result_Fitness = 10 - result_Fitness;  

 %%%%%%%%%%%%%%%更新全局变量%%%%%%%%%%%%%%%%%%%%
 fitnesstimes = fitnesstimes + 1;
 if fitnesstimes == (size_particle + 1)
    fitnesstimes = 1;
 end 
 if fitnesstimes == 1                 
    [Ytemp,MatchGenesWithParticles] = sort(rand(1,size_particle)); 
 end 