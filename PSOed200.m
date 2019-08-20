%Author: Luo Qiang
%Abstract: minimizing function by PSO
%version: 1.0.1
%date: 7 June, 2007
%���룺 
%   Ŀ�꺯�������� ObjPara = struct('NCBP',[min, max], 'Dim', dim_objfunction,'Name', objfunctionname)
%   ��ѡ������ 'Schaffer','Rastrigian','Rossenbrock','Six-humpCamelback',
%   'Sphere','2n-minima','Griewank'.��myFun.m�ж��塣
%    �㷨���������� AlgPara = struct('MaxIte',max_iteration, 'WeightType', flag_weight, 'MatchType', flag_gamatch,
%                               'IniRange',[ini_low, ini_high], 'Algorithm', flag_algorithm, 'Size', size_particle);
%       �Ѿ�ʵ�ֵ��㷨�У� SPSO, LPSO, FPSO, RPSO
%       ����������max_iteration,
%       �������Ӳ��õ�һȨ��:flag_weight = 1 ����ÿ�����Ӳ�ͬ��Ȩ��flag_weight = 2; 
%       flag_weightcontrol  ����ʲô�������ƹ���Ȩ��: LPSO,FPSO,RPSO
%       ����RPSO��Ҫ���ǻ��������ӵ�ƥ�䷽���������Ҫ�������ƥ��: flag_gamatch = 1; 
%           �ʶȶ�ƥ�� = 0; ����ʽѡ��ƥ�� = 2��   
%       ��ʼ����Χ�� ini_low_particle, ini_high_particle 
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
errorBand = 1e-4; %����4λ��Ч����


%ÿһά�ϵ�X����Χ
abs_max_particle = 10;
%ÿһά�ϵ��ٶȵ�����Χ

abs_max_V = 10;
%ÿһά�ϵ��ٶȵĳ�ʼ����Χ
ini_low_V = 0;
ini_high_V = 10;
%Ȩ�صĳ�ʼ������
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
    %ÿ�����ӣ���ͬ������ÿ�����ӵ�ÿһά����ͬ��
    ini_inertia_weight_high = 0.90;
    ini_inertia_weight_low = 0.40;
    W(:,:,:) = rand(size_particle,1) * ones(1,dim_objfunction)...
        * (ini_inertia_weight_high ...
        - ini_inertia_weight_low) + ini_inertia_weight_low; 
    ini_W = W(:,:,1);
end

%Compute the ObjFunction Value for each Particle
F(:,iteration) = ComputeObjFunctionForParticle(Particles(:,:,iteration), name_objfunction);

%ÿ�����ӵľֲ���Сֵ
LBest(:,:) = F(:,iteration);
%��¼�ֲ���Сλ��
local_minimum_particle(:,:,:) =  Particles(:,:,iteration);

%�������ӵ�ȫ����Сֵ
[GBest(iteration),GIndex(iteration)] = min(LBest(:,iteration));  
%��¼ȫ����Сλ��Ϊ�˸��·��㣬�����γ�һ��ÿһҳ��ÿһ�ж���ͬ����ά����
global_minimum_particle(:,:,:) = ones(size_particle,1) * Particles(GIndex(iteration),:,iteration);

%%��¼GA���ݻ�����
GACourse(iteration) = 0;
%%��¼GA����Ӧ��
GAFitness(:,iteration) = zeros(1,size_particle);
%%��¼�ڵ�һȨ���㷨�У�ÿ��ѡ�еĹ�������
CG_Index = 0;
%%��¼��fire�������Ӧ�Ⱥ���ֵ
Index_Fitness = 0;
%-----------------------------------------%



%-----------------------------------------%






%%��ͼ
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

%��������
no_decrease_count = 0;
term_condition = 1;
while(term_condition)    

    %update velocity
    V(:,:,iteration+1) = W(:,:,iteration) .* V(:,:,iteration)...
        - rand * (c1_local(:,iteration) * ones(1,dim_objfunction)) ...
        .* ( Particles(:,:,iteration) - local_minimum_particle(:,:,iteration) ) ...
        - rand * (c2_global(:,iteration) * ones(1,dim_objfunction))...
        .* ( Particles(:,:,iteration) - global_minimum_particle(:,:,iteration) );     
    %Խ����V
    V(:,:,iteration+1) = min(V(:,:,iteration+1),abs_max_V);  %��ȫ���ٶ�С��abs_max_V
    V(:,:,iteration+1) = max(V(:,:,iteration+1),-abs_max_V); %��ȫ���ٶȴ���-abs_max_V
    
    %particle fly
    Particles(:,:,iteration+1) = Particles(:,:,iteration) + V(:,:,iteration+1);
    %Խ����X
    Particles(:,:,iteration+1) = min(Particles(:,:,iteration+1),abs_max_particle);
    Particles(:,:,iteration+1) = max(Particles(:,:,iteration+1),-abs_max_particle);
    
%     %ÿ�ε���������һ��������������ӣ��Ա��ֶ�����
%     Particles(1:5,:,iteration+1) = ini_low_particle ...
%     + (ini_high_particle - ini_low_particle) ...
%     * rand(5, dim_objfunction, 1);
    
    
    %Compute the ObjFunction Value for each Particle
    F(:,iteration+1) = ComputeObjFunctionForParticle(Particles(:,:,iteration+1), name_objfunction);
  
    %update local best 
    for i = 1:size_particle
        if ( F(i,iteration+1) < LBest(i,iteration) )
            %�������ֵС�ˣ�����            
            LBest(i,iteration+1) = F(i,iteration+1);
            local_minimum_particle(i,:,iteration+1) = Particles(i,:,iteration+1);            
        else
            %�������ֵû�б�С��������            
            LBest(i,iteration+1) = LBest(i,iteration);
            local_minimum_particle(i,:,iteration+1) = local_minimum_particle(i,:,iteration);            
        end
    end
    
    %update global best 
    [GBest(iteration+1),GIndex(iteration+1)] = min(LBest(:,iteration+1));  
    if ( GBest(iteration+1) < GBest(iteration) )
        %�����С������
        global_minimum_particle(:,:,iteration+1) = ones(size_particle,1) * Particles(GIndex(iteration+1),:,iteration+1);        
    else
        %���򣬱���ԭ����ȫ����С
        global_minimum_particle(:,:,iteration+1) = global_minimum_particle(:,:,iteration);
        GBest(iteration+1) = GBest(iteration);
    end
    %��¼��С��������
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
            %% �㶨Ȩ��
            %Method Ia��SimplePSO with uniW
            W(:,:,iteration+1) = W(:,:,iteration);
            %Method Ib: SimplePSO with divW
            %%---------------------------------------------------------------------           
        case 'LPSO'
             %%%%%%%
     %%%%%%%
     plotkey = 0;
     %%%%%%%
     %%%%%%%
            %%���Եݼ�Ȩ��   
            if (flag_weight == 1)
                %Method IIa��LinearPSO
                W(:,:,iteration+1) = uniUpdateWeight(iteration+1,max_iteration, size_particle,dim_objfunction); %LPSO     
            else
                %Method IIb��LinearPSO
                W(:,:,iteration+1) = divUpdateWeight(iteration+1,max_iteration, ini_W);
            end           
             %%------------------------------------------------------------
        case 'FPSO'
             %%%%%%%
     %%%%%%%
     plotkey = 0;
     %%%%%%%
     %%%%%%%
            %%ģ������Ȩ��           
            if (flag_weight == 1)
                %Method IIIa: FuzzyPSO
                CBPE = [NCBP_min, NCBP_max, GBest(iteration+1)];   %Ԥ�ں�����Сֵ�����������ֵ����ǰ������ֵ
                W_next = uniFAUpdateWeight(CBPE, W(1,1,iteration)); %����Ȩ�ؾ���Ӧ�þ���һ����
                W(:,:,iteration+1) = W_next * ones(size_particle, dim_objfunction, 1);
            else
                %Method IIIb:
                CBPE = [NCBP_min, NCBP_max, GBest(iteration+1)]; 
                W(:,:,iteration+1) = divFAUpdateWeight(CBPE, W(:,:,iteration));%����Ȩ�ؾ���ÿһ�ж���ͬ����ÿ�����ӵ�ÿһά����һ��Ȩ��
            end            
            %%-----------------------------------------------------------
            
        case 'RPSO'
            %%%%%%%
            %%%%%%%
            plotkey = 1;
            %%%%%%%
            %%%%%%%

            %%ģ���Ļ��㷨����Ȩ��
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
                %������Ӧ�Ⱥ���
                if iteration == 1
                    %�����ʼCG����Ӧ��
                    %��ÿ������fireһ�Σ��ɴ�����������Ӧ��
                    GAFitness(:,1) = InitilizeFitness(iniPop,Fit_paras);
                else
                    %����fire�������ֵ�ı��˶���
                    GAFitness(:,iteration) = LBest(:,iteration) - F(:,iteration+1);                                      
                end
                %���滯��[-min,max],����������Ӧ�ȶ���      
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
                   GAFitness(:,iteration) =  exp(GAFitness(:,iteration)); %����Griewank�����ĺ���������ֵ��С

               
                                
                %���ݵ�ǰ����Ӧ�Ⱥ�������������ȺnewPop             
                if (flag_gamatch == 1)                    
                    %���������CGִ���Ŵ�����,indexMatrix�ĵڼ��ж�Ӧ�˵ڼ�������ƥ���������������                    
                    [Ytemp, indexMatrix] = sort(rand(size_particle,2));
                    %�����µ���Ⱥ
                    NewPop = GeneticOperation(OldPop, indexMatrix);
                elseif (flag_gamatch == 0)
                    %�̶�����������CGִ���Ŵ�����
                    indexMatrix = [1:size_particle;2:size_particle,1]';
                    %�����µ���Ⱥ
                    NewPop = GeneticOperation(OldPop, indexMatrix);                    
                else               
                    %����ѡ������CGִ���Ŵ�����
                    indexMatrix = RouletteSelect(GAFitness(:,iteration));                    
                    %�����µ���Ⱥ
                    NewPop = GeneticOperation(OldPop, indexMatrix);                
                end

                CBPE = [NCBP_min, NCBP_max, GBest(iteration+1)];
                %����NewPop���¹���Ȩ��,�Ѿ������ƥ��
                W(:,:,iteration+1) = divFAUpdateWeight_RPSO(CBPE, W(:,:,iteration), NewPop); 
            end
            %%-----------------------------------------------------------
            %%����OldPop
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
             
    %%%������ͼ����������������������������������������
    %%%%%
    %�������ŵķ�������ʵ�˶�����ļн�,����ֵ��0
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
%     %�������ӷ��е����   
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
    % ¼��ʹ��
%     if (iteration == 2) 
%         pause(10); 
%     end
    %%%����������������������������������������������������   
        
end

%%%-----��ͼ,��ƫ�ǽ���ͳ�Ʒ���,ƽ��ƫ�ǵı仯
% figure('Name','Angle')
% plot(1:size(theta,2), mean(theta));
%%%-------------------------------------------

Find = struct('best',GBest(iteration-1),...
    'BestGA',max(GAFitness), 'AveGA', mean(GAFitness),...
    'BestPSO', GBest(1:iteration-1), 'AvePSO', mean(LBest(:,1:iteration-1))); 

