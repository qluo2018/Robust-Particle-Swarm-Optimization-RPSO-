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
%      'sLBest_Curr', LBest(:,iteration),'sW_Current',W(1,1,iteration))
%�����result_Fitness���Ļ��������Ӧ��
%%%% %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%% 
%%%% ÿ������������1��1���ƥ��
%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%% 

function result_Fitness = unimyFitness_Match ( x, Fit_paras )

%%%%%%%%%%%%%%%% ����ȫ�ֱ�������ƥ����������
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


%%%%%%%%%%%%%%%%ʹ��ȫ�ֱ���%%%%%%%%%%%%%%%%%%%%%%%
index_temp = MatchGenesWithParticles(fitnesstimes); 


%update W
 CBPE = [Fit_paras.sNCBP(1), Fit_paras.sNCBP(2), Fit_paras.sGBest_Curr];
 W_NextIter = uniFAUpdateWeight_FCAPSO(CBPE, Fit_paras.sW_Current, delta_g); %����Ȩ�ؾ���Ӧ�þ���һ����
 
 clear V_Next;
 clear Particles_Next;
 %update velocity
 V_Next = W_NextIter * Fit_paras.sV_Curr(index_temp,:)...
        - rand * Fit_paras.sc1_local_current(1) ...
        * ( Fit_paras.sParticles_Curr(index_temp,:) - Fit_paras.slocal_minimum_particle_Curr(index_temp,:) ) ...
        - rand * Fit_paras.sc2_global_current(1)...
        * ( Fit_paras.sParticles_Curr(index_temp,:) - Fit_paras.sglobal_minimum_particle_Curr(index_temp,:) );     

    %Խ����V
    V_Next = min(V_Next,Fit_paras.sabs_max_V);  %��ȫ���ٶ�С��abs_max_V
    V_Next = max(V_Next,-Fit_paras.sabs_max_V); %��ȫ���ٶȴ���-abs_max_V
    
    %particle fly
    Particles_Next =  Fit_paras.sParticles_Curr (index_temp,:) + V_Next;
    %Խ����X
    Particles_Next = min(Particles_Next, Fit_paras.sabs_max_particle);
    Particles_Next = max(Particles_Next, -Fit_paras.sabs_max_particle);
    

    
    %Compute the ObjFunction Value for each Particle
    F_Next = ComputeObjFunctionForParticle(Particles_Next);
    

   
       
   result_Fitness = 0;
    %�þֲ�����ֵ����ֵ��ֲ����űȽϣ�������Ӧ��
        if ( F_Next < Fit_paras.sLBest_Curr(index_temp) )
            %�������ֵС�ˣ�����            
            result_Fitness = 1 +  W_NextIter;                              
        end
    %���������űȽϣ�������Ӧ��      
    if ( F_Next < Fit_paras.sGBest_Curr )
        %�����С������
        result_Fitness = 3;       
    end
    result_Fitness = 10 - result_Fitness;  

 %%%%%%%%%%%%%%%����ȫ�ֱ���%%%%%%%%%%%%%%%%%%%%
 fitnesstimes = fitnesstimes + 1;
 if fitnesstimes == (size_particle + 1)
    fitnesstimes = 1;
 end 
 if fitnesstimes == 1                 
    [Ytemp,MatchGenesWithParticles] = sort(rand(1,size_particle)); 
 end 