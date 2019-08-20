%%�Գ�ʼ��Ⱥ������Ӧ��
%%ֻ����ÿ������fireһ��
%%���룺��ʼ��Ⱥ����Ӧ�Ȳ���
%%�������ʼ��Ⱥ����Ӧ��


function IniFitness = InitilizeFitness(iniPop,Fit_paras)

size_particle = size(iniPop);
dim_objfunction = Fit_paras.sSizePSO(2);


%ÿ��delta_g��Ӧ��һ��CG,CG�����ӵĶ�Ӧ��ϵ�Ѿ���Ӧ����
%update W
CBPE = [Fit_paras.sNCBP(1), Fit_paras.sNCBP(2), Fit_paras.sGBest_Curr]; %����NCBPEʱ���ǻ������������Ⱥ�����㣬���ټ�������ֻ��ÿ�����ӵĹ���Ȩ�ز�һ��
W_NextIter = divFAUpdateWeight_RPSO(CBPE, Fit_paras.sW_Current, iniPop); %����Ȩ�ؾ���ÿһ�ж���ͬ
clear V_Next;
clear Particles_Next;

 for index_temp = 1: size_particle
 %update velocity
 V_Next = W_NextIter(index_temp,1) * Fit_paras.sV_Curr(index_temp,:)...
        - rand * Fit_paras.sc1_local_current(1) ...
        * ( Fit_paras.sParticles_Curr(index_temp,:) - Fit_paras.slocal_minimum_particle_Curr(index_temp,:) ) ...
        - rand * Fit_paras.sc2_global_current(1)...
        * ( Fit_paras.sParticles_Curr(index_temp,:) - Fit_paras.sglobal_minimum_particle_Curr(index_temp,:) );     
 %Խ����V
 V_Next = min(V_Next,Fit_paras.sabs_max_V);  %��ȫ���ٶ�С��abs_max_V
 V_Next = max(V_Next,-Fit_paras.sabs_max_V); %��ȫ���ٶȴ���-abs_max_V
    
    
 %particle fly
 Particles_Next =  Fit_paras.sParticles_Curr(index_temp,:) + V_Next;
 %Խ����X
 Particles_Next = min(Particles_Next, Fit_paras.sabs_max_particle);
 Particles_Next = max(Particles_Next, -Fit_paras.sabs_max_particle);
 
 %Compute the ObjFunction Value for each Particle
 F_Next = ComputeObjFunctionForParticle(Particles_Next,Fit_paras.sName);
 
 %�����ʼ��Ӧ��
 %����Ŀ�꺯��ֵ�����˶��٣�
 IniFitness(index_temp) = Fit_paras.sLBest_Curr(index_temp) - F_Next;
 
 end

