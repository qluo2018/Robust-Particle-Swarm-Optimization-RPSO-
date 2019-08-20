%%ͳ�Ʒ���PSO������
%%ÿ���㷨����30�Σ������������ݣ��Լ�ƽ��ֵ�ͱ�׼����Ҷ�ƽ��ֵ��ͼ
clear
clc
% %�㶨ͳһȨ��
% 
max_i = 4;
max_repeat = 30;
max_iteration = 300;
function_name = 'Schaffer';

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


%��ʱֻ���������Ժ���
%Rastrigian [0,70]
%Rossenbrock [0,500]
%Griewank [0,0.05]
%Sphere [0,infinite]
%Schaffer [0,infinite]
%�Ҷ�����0��100




clear ObjPara;
ObjPara = struct('NCBP',[0, 10], 'Dim', 2, 'Name', function_name);
clear AlgPara;   
AlgPara = struct('MaxIte',max_iteration, 'WeightType', 1, ...
            'MatchType', 1, 'IniRange',[-5, +5],...
            'Algorithm', 'LPSO', 'Size', 20);
%%LPSO
for i = 1:max_i
    dim_obj = 20 + 60 * (i-1);
    ObjPara.Dim = dim_obj;    
    for repeat_times = 1 : max_repeat       
        % ÿһ�б�ʾһ���߶�
        % ÿһ�б�ʾһ���ظ�
        AlgPara.WeightType = 1; %����Ϊÿ������ʹ����ͬ��Ȩ��
        t = cputime;        
        Find_uni_LPSO = PSOed200(ObjPara,AlgPara);     
        t_uni_LPSO(repeat_times,i) = cputime- t;
        Best_uni_LPSO(repeat_times,i) = Find_uni_LPSO.best; 
        BestPSO_uni_LPSO(:,repeat_times,i) = Find_uni_LPSO.BestPSO;
        AvePSO_uni_LPSO(:,repeat_times,i) = Find_uni_LPSO.AvePSO;
                
        AlgPara.WeightType = 2; %����Ϊÿ������ʹ�ò�ͬ��Ȩ��
        t = cputime;        
        Find_div_LPSO = PSOed200(ObjPara,AlgPara);     
        t_div_LPSO(repeat_times,i) = cputime- t;
        Best_div_LPSO(repeat_times,i) = Find_div_LPSO.best; 
        BestPSO_div_LPSO(:,repeat_times,i) = Find_div_LPSO.BestPSO;
        AvePSO_div_LPSO(:,repeat_times,i) = Find_div_LPSO.AvePSO;       
        
    end
end

%%д�ļ�����ֵ
fileWrite('uniLPSO.txt',Best_uni_LPSO, max_repeat, max_i);
fileWrite('divLPSO.txt',Best_div_LPSO, max_repeat, max_i);

%%ƽ���������̻�ͼ
for i = 1:max_i    
    plotAveCourse(BestPSO_uni_LPSO(:,:,i),AvePSO_uni_LPSO(:,:,i),...
        BestPSO_div_LPSO(:,:,i), AvePSO_div_LPSO(:,:,i),...
        'LPSO','uniweight LPSO Course', 'diversity LPSO Course',i); 
end



%%RPSO
clear AlgPara.Algorithm;
AlgPara.Algorithm = 'RPSO';
AlgPara.WeightType = 2; %����Ϊÿ������ʹ�ò�ͬ��Ȩ��
for i = 1:max_i
    dim_obj = 20 + 60 * (i-1);
    ObjPara.Dim = dim_obj;     
    for repeat_times = 1 : max_repeat
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        AlgPara.MatchType = 0;% �̶�ƥ��
        clear Find_fixed_RPSO;
        t = cputime;
        Find_fixed_RPSO = PSOed200(ObjPara,AlgPara);  
        t_fixed_RPSO(repeat_times,i) = cputime - t;
        
        Best_fixed_RPSO(repeat_times,i) = Find_fixed_RPSO.best;
        BestGA_fixed_RPSO(:,repeat_times,i) = Find_fixed_RPSO.BestGA;
        AveGA_fixed_RPSO(:,repeat_times,i) = Find_fixed_RPSO.AveGA;
        BestPSO_fixed_RPSO(:,repeat_times,i) = Find_fixed_RPSO.BestPSO;
        AvePSO_fixed_RPSO(:,repeat_times,i) = Find_fixed_RPSO.AvePSO;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        AlgPara.MatchType = 1;  %���ƥ��
        clear Find_random_RPSO;
        t = cputime;
        Find_random_RPSO = PSOed200(ObjPara,AlgPara);
        t_random_RPSO(repeat_times,i) = cputime - t; 
        
        Best_random_RPSO(repeat_times,i) = Find_random_RPSO.best;
        BestGA_random_RPSO(:,repeat_times,i) = Find_random_RPSO.BestGA;
        AveGA_random_RPSO(:,repeat_times,i) = Find_random_RPSO.AveGA;
        BestPSO_random_RPSO(:,repeat_times,i) = Find_random_RPSO.BestPSO;
        AvePSO_random_RPSO(:,repeat_times,i) = Find_random_RPSO.AvePSO; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        AlgPara.MatchType = 2; %����Ӧ������ѡ��
        clear Find_roulette_RPSO;
        t = cputime;
        Find_roulette_RPSO = PSOed200(ObjPara,AlgPara);
        t_roulette_RPSO(repeat_times,i) = cputime - t; 
        
        Best_roulette_RPSO(repeat_times,i) = Find_roulette_RPSO.best;        
        BestGA_roulette_RPSO(:,repeat_times,i) = Find_roulette_RPSO.BestGA;
        AveGA_roulette_RPSO(:,repeat_times,i) = Find_roulette_RPSO.AveGA;
        BestPSO_roulette_RPSO(:,repeat_times,i) = Find_roulette_RPSO.BestPSO;
        AvePSO_roulette_RPSO(:,repeat_times,i) = Find_roulette_RPSO.AvePSO; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
end

%%д�ļ�����ֵ
fileWrite('fixedRPSO.txt',Best_fixed_RPSO, max_repeat, max_i);
fileWrite('randomRPSO.txt',Best_random_RPSO, max_repeat, max_i);
fileWrite('rouletteRPSO.txt',Best_roulette_RPSO, max_repeat, max_i);

%%ƽ���������̻�ͼ
for i = 1:max_i    
    plotAveCourse(BestGA_fixed_RPSO(:,:,i),AveGA_fixed_RPSO(:,:,i),...
        BestPSO_fixed_RPSO(:,:,i), AvePSO_fixed_RPSO(:,:,i),...
        'fixed_RPSO','GA Course', 'PSO Course',i); 
    plotAveCourse(BestGA_random_RPSO(:,:,i),AveGA_random_RPSO(:,:,i),...
        BestPSO_random_RPSO(:,:,i), AvePSO_random_RPSO(:,:,i),...
        'random_RPSO','GA Course', 'PSO Course',i);  
    plotAveCourse(BestGA_roulette_RPSO(:,:,i),AveGA_roulette_RPSO(:,:,i),...
        BestPSO_roulette_RPSO(:,:,i), AvePSO_roulette_RPSO(:,:,i),...
        'roulette_RPSO','GA Course', 'PSO Course',i);    
end


clear AlgPara.Algorithm;
AlgPara.Algorithm = 'FPSO';
%%FPSO
for i = 1:max_i
    dim_obj = 20 + 60 * (i-1);
    ObjPara.Dim = dim_obj;  
    for repeat_times = 1 : max_repeat
        % ÿһ�б�ʾһ���߶�
        % ÿһ�б�ʾһ���ظ�
        AlgPara.WeightType = 1; %����Ϊÿ������ʹ����ͬ��Ȩ��
        t = cputime;
        Find_uni_FPSO = PSOed200(ObjPara,AlgPara);     
        t_uni_FPSO(repeat_times,i) = cputime- t;
        Best_uni_FPSO(repeat_times,i) = Find_uni_FPSO.best; 
        BestPSO_uni_FPSO(:,repeat_times,i) = Find_uni_FPSO.BestPSO;
        AvePSO_uni_FPSO(:,repeat_times,i) = Find_uni_FPSO.AvePSO;
        
        AlgPara.WeightType = 2; %����Ϊÿ������ʹ�ò�ͬ��Ȩ��
        t = cputime;        
        Find_div_FPSO = PSOed200(ObjPara,AlgPara);     
        t_div_FPSO(repeat_times,i) = cputime- t;
        Best_div_FPSO(repeat_times,i) = Find_div_FPSO.best; 
        BestPSO_div_FPSO(:,repeat_times,i) = Find_div_FPSO.BestPSO;
        AvePSO_div_FPSO(:,repeat_times,i) = Find_div_FPSO.AvePSO;   
    end
end
%%д�ļ�����ֵ
fileWrite('uniFPSO.txt',Best_uni_FPSO, max_repeat, max_i);
fileWrite('divFPSO.txt',Best_div_FPSO, max_repeat, max_i);

%%ƽ���������̻�ͼ
for i = 1:max_i    
    plotAveCourse(BestPSO_uni_FPSO(:,:,i),AvePSO_uni_FPSO(:,:,i),...
        BestPSO_div_FPSO(:,:,i), AvePSO_div_FPSO(:,:,i),...
        'FPSO','uniweight FPSO Course', 'diversity FPSO Course',i); 
end

%%д�ļ�ƽ������ʱ��
fid = fopen('ProcessTime.txt','wt');
for i = 1 : max_i    
    Ave_t_uni_LPSO  = mean(t_uni_LPSO);
    Ave_t_div_LPSO  = mean(t_div_LPSO);
    Ave_t_uni_FPSO  = mean(t_uni_FPSO);
    Ave_t_div_FPSO  = mean(t_div_FPSO);    
    Ave_t_fixed_RPSO  = mean(t_fixed_RPSO);
    Ave_t_random_RPSO  = mean(t_random_RPSO);
    Ave_t_roulette_RPSO  = mean(t_roulette_RPSO);    
    fprintf(fid, '%12.4f  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f\n',...
       Ave_t_uni_LPSO(i), Ave_t_div_LPSO(i),Ave_t_uni_FPSO(i),...
       Ave_t_div_FPSO(i), Ave_t_fixed_RPSO(i), Ave_t_random_RPSO(i),...
       Ave_t_roulette_RPSO(i));
end
fclose(fid);