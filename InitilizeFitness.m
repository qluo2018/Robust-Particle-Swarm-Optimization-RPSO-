%%对初始种群计算适应度
%%只能让每条规则都fire一次
%%输入：初始种群和适应度参数
%%输出：初始种群的适应度


function IniFitness = InitilizeFitness(iniPop,Fit_paras)

size_particle = size(iniPop);
dim_objfunction = Fit_paras.sSizePSO(2);


%每个delta_g对应了一个CG,CG与粒子的对应关系已经对应好了
%update W
CBPE = [Fit_paras.sNCBP(1), Fit_paras.sNCBP(2), Fit_paras.sGBest_Curr]; %计算NCBPE时我们还是针对整个种群来计算，减少计算量，只是每个粒子的惯性权重不一样
W_NextIter = divFAUpdateWeight_RPSO(CBPE, Fit_paras.sW_Current, iniPop); %整个权重矩阵每一列都相同
clear V_Next;
clear Particles_Next;

 for index_temp = 1: size_particle
 %update velocity
 V_Next = W_NextIter(index_temp,1) * Fit_paras.sV_Curr(index_temp,:)...
        - rand * Fit_paras.sc1_local_current(1) ...
        * ( Fit_paras.sParticles_Curr(index_temp,:) - Fit_paras.slocal_minimum_particle_Curr(index_temp,:) ) ...
        - rand * Fit_paras.sc2_global_current(1)...
        * ( Fit_paras.sParticles_Curr(index_temp,:) - Fit_paras.sglobal_minimum_particle_Curr(index_temp,:) );     
 %越界检查V
 V_Next = min(V_Next,Fit_paras.sabs_max_V);  %让全部速度小于abs_max_V
 V_Next = max(V_Next,-Fit_paras.sabs_max_V); %让全部速度大于-abs_max_V
    
    
 %particle fly
 Particles_Next =  Fit_paras.sParticles_Curr(index_temp,:) + V_Next;
 %越界检查X
 Particles_Next = min(Particles_Next, Fit_paras.sabs_max_particle);
 Particles_Next = max(Particles_Next, -Fit_paras.sabs_max_particle);
 
 %Compute the ObjFunction Value for each Particle
 F_Next = ComputeObjFunctionForParticle(Particles_Next,Fit_paras.sName);
 
 %计算初始适应度
 %计算目标函数值降低了多少，
 IniFitness(index_temp) = Fit_paras.sLBest_Curr(index_temp) - F_Next;
 
 end

