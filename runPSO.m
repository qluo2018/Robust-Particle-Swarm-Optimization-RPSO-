%%统计分析PSO的性能
%%每种算法运行30次，保留运行数据，以及平均值和标准差，并且对平均值画图
clear
clc
% %恒定统一权重
% 
max_i = 4;
max_repeat = 30;
max_iteration = 300;
function_name = 'Schaffer';

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


%暂时只用三个测试函数
%Rastrigian [0,70]
%Rossenbrock [0,500]
%Griewank [0,0.05]
%Sphere [0,infinite]
%Schaffer [0,infinite]
%我都采用0，100




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
        % 每一列表示一个尺度
        % 每一行表示一次重复
        AlgPara.WeightType = 1; %设置为每个粒子使用相同的权重
        t = cputime;        
        Find_uni_LPSO = PSOed200(ObjPara,AlgPara);     
        t_uni_LPSO(repeat_times,i) = cputime- t;
        Best_uni_LPSO(repeat_times,i) = Find_uni_LPSO.best; 
        BestPSO_uni_LPSO(:,repeat_times,i) = Find_uni_LPSO.BestPSO;
        AvePSO_uni_LPSO(:,repeat_times,i) = Find_uni_LPSO.AvePSO;
                
        AlgPara.WeightType = 2; %设置为每个粒子使用不同的权重
        t = cputime;        
        Find_div_LPSO = PSOed200(ObjPara,AlgPara);     
        t_div_LPSO(repeat_times,i) = cputime- t;
        Best_div_LPSO(repeat_times,i) = Find_div_LPSO.best; 
        BestPSO_div_LPSO(:,repeat_times,i) = Find_div_LPSO.BestPSO;
        AvePSO_div_LPSO(:,repeat_times,i) = Find_div_LPSO.AvePSO;       
        
    end
end

%%写文件最优值
fileWrite('uniLPSO.txt',Best_uni_LPSO, max_repeat, max_i);
fileWrite('divLPSO.txt',Best_div_LPSO, max_repeat, max_i);

%%平均迭代过程画图
for i = 1:max_i    
    plotAveCourse(BestPSO_uni_LPSO(:,:,i),AvePSO_uni_LPSO(:,:,i),...
        BestPSO_div_LPSO(:,:,i), AvePSO_div_LPSO(:,:,i),...
        'LPSO','uniweight LPSO Course', 'diversity LPSO Course',i); 
end



%%RPSO
clear AlgPara.Algorithm;
AlgPara.Algorithm = 'RPSO';
AlgPara.WeightType = 2; %设置为每个粒子使用不同的权重
for i = 1:max_i
    dim_obj = 20 + 60 * (i-1);
    ObjPara.Dim = dim_obj;     
    for repeat_times = 1 : max_repeat
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        AlgPara.MatchType = 0;% 固定匹配
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
        AlgPara.MatchType = 1;  %随机匹配
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
        AlgPara.MatchType = 2; %按适应度轮盘选择
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

%%写文件最优值
fileWrite('fixedRPSO.txt',Best_fixed_RPSO, max_repeat, max_i);
fileWrite('randomRPSO.txt',Best_random_RPSO, max_repeat, max_i);
fileWrite('rouletteRPSO.txt',Best_roulette_RPSO, max_repeat, max_i);

%%平均迭代过程画图
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
        % 每一列表示一个尺度
        % 每一行表示一次重复
        AlgPara.WeightType = 1; %设置为每个粒子使用相同的权重
        t = cputime;
        Find_uni_FPSO = PSOed200(ObjPara,AlgPara);     
        t_uni_FPSO(repeat_times,i) = cputime- t;
        Best_uni_FPSO(repeat_times,i) = Find_uni_FPSO.best; 
        BestPSO_uni_FPSO(:,repeat_times,i) = Find_uni_FPSO.BestPSO;
        AvePSO_uni_FPSO(:,repeat_times,i) = Find_uni_FPSO.AvePSO;
        
        AlgPara.WeightType = 2; %设置为每个粒子使用不同的权重
        t = cputime;        
        Find_div_FPSO = PSOed200(ObjPara,AlgPara);     
        t_div_FPSO(repeat_times,i) = cputime- t;
        Best_div_FPSO(repeat_times,i) = Find_div_FPSO.best; 
        BestPSO_div_FPSO(:,repeat_times,i) = Find_div_FPSO.BestPSO;
        AvePSO_div_FPSO(:,repeat_times,i) = Find_div_FPSO.AvePSO;   
    end
end
%%写文件最优值
fileWrite('uniFPSO.txt',Best_uni_FPSO, max_repeat, max_i);
fileWrite('divFPSO.txt',Best_div_FPSO, max_repeat, max_i);

%%平均迭代过程画图
for i = 1:max_i    
    plotAveCourse(BestPSO_uni_FPSO(:,:,i),AvePSO_uni_FPSO(:,:,i),...
        BestPSO_div_FPSO(:,:,i), AvePSO_div_FPSO(:,:,i),...
        'FPSO','uniweight FPSO Course', 'diversity FPSO Course',i); 
end

%%写文件平均运行时间
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