%��Ӧ�Ⱥ���
%����������
%���룺 
% GA ����С��

%Fit_paras =
%struct('sNCBP',[NCBP_min,NCBP_max],'sGBest_Curr',GBest(iteration),'sSizePSO',[size_particles,dim_objfunction],...
%      'sc1_local_current',c1_local_current(:,iteration),'sc2_global_current',c2_global_current(:,iteration), ...
%      'sV_Curr', V_Current(:,:,iteration), 'sParticles_Curr', Particles(:,:,iteration), ...
%      'slocal_minimum_particle_Curr', local_minimum_particle(:,:,iteration) , ...
%      'sglobal_minimum_particle_Curr', global_minimum_particle(:,:,iteration), ... 
%      'sabs_max_V', abs_max_V,'sabs_max_particle',abs_max_particle,...
%      'sLBest_Curr', LBest(:,iteration),'sW_Current',W(:,:,iteration))
%�����result_Fitness���Ļ��������Ӧ��
%%%% %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%% 
%%%% ÿ���������Ӧ���Ƕ���������Ⱥ���Ч�������ѡ����õĻ�������Ϊ��һ��ʹ�õĹ������
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
 CBPE = [Fit_paras.sNCBP(1), Fit_paras.sNCBP(2), Fit_paras.sGBest_Curr]; %����NCBPEʱ���ǻ������������Ⱥ�����㣬���ټ�������ֻ��ÿ�����ӵĹ���Ȩ�ز�һ��
 W_NextIter = divFAUpdateWeight_FCAPSO(CBPE, Fit_paras.sW_Current, delta_g); %����Ȩ�ؾ���ÿһ�ж���ͬ
 
%update velocity
    V_Next(:,:) = W_NextIter .* Fit_paras.sV_Curr...
        - rand * (Fit_paras.sc1_local_current * ones(1,dim_objfunction)) ...
        .* ( Fit_paras.sParticles_Curr - Fit_paras.slocal_minimum_particle_Curr ) ...
        - rand * (Fit_paras.sc2_global_current * ones(1,dim_objfunction))...
        .* ( Fit_paras.sParticles_Curr - Fit_paras.sglobal_minimum_particle_Curr );     
    %Խ����V
    V_Next = min(V_Next,Fit_paras.sabs_max_V);  %��ȫ���ٶ�С��abs_max_V
    V_Next = max(V_Next,-Fit_paras.sabs_max_V); %��ȫ���ٶȴ���-abs_max_V
    
    %particle fly
    Particles_Next(:,:) =  Fit_paras.sParticles_Curr + V_Next;
    %Խ����X
    Particles_Next = min(Particles_Next, Fit_paras.sabs_max_particle);
    Particles_Next = max(Particles_Next, -Fit_paras.sabs_max_particle);
    
    
    
    %Compute the ObjFunction Value for each Particle
    F_Next = ComputeObjFunctionForParticle(Particles_Next);
    
    %�þֲ�����ֵ����ֵ��ֲ����űȽϣ�������Ӧ��
    result_Fitness = 0;
    for i = 1:size_particle
        if ( F_Next(i) < Fit_paras.sLBest_Curr(i) )
            %�������ֵС�ˣ�����            
            result_Fitness = result_Fitness + 1 + W_NextIter(i,1);   
             %�������ֵС�ˣ�����            
            LBest_Next(i) = F_Next(i);             
        else
            %�������ֵû�б�С��������            
            %result_Fitness = result_Fitness + 0;    
             %�������ֵû�б�С��������            
            LBest_Next(i) = Fit_paras.sLBest_Curr(i);            
        end
    end
    result_Fitness = result_Fitness / size_particle;   %�Ļ����������Ⱥ��Ӧ���޸ĵ�ƽ��Ч��
    
    %���ȫ������ֵ��ø��£���ô����������Ӧ�ȱȽϸ� 
    [GBest_Next,GIndex_Next] = min(abs(LBest_Next));  
    if ( GBest_Next < Fit_paras.sGBest_Curr )
        %�����С������
        result_Fitness = 3;       
    end
   result_Fitness = 10 - result_Fitness;  
 
 

