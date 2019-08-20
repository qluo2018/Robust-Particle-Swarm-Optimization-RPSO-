%%统计分析PSO的性能
%%每种算法运行30次，保留运行数据，以及平均值和标准差，并且对平均值画图
clear
clc
% %恒定统一权重
% 
global sigma_noise;
max_i = 6;
max_repeat = 30;
max_iteration = 300;
%%%max_iteration = 1000;
dim_obj = 2; %50,100

%%%FCAPSO_match
% for i = 1:max_i
%     sigma_noise = 2*(i-1)/10;
%     for repeat_times = 1 : max_repeat
%         %t = cputime;
%         %Best_FCAPSO_match_uniweight(repeat_times,i)= PSOed100(dim_obj,max_iteration,1,'FCAPSO',1);   %1 match
%         %t_FCAPSOmu(i) = cputime - t;
%         %t = cputime;
%         Best_FCAPSO_match_divweight(repeat_times,i) = PSOed100(dim_obj,max_iteration,2,'FCAPSO',1);
%         %t_FCAPSOmd(i) = cputime - t;    
%     end
% end
% %%写文件
% % fid = fopen('n_uFCAPSO_match_30.txt','wt');
% % for repeat_times = 1:max_repeat   
% %     for i = 1:max_i
% %         fprintf(fid,'%12.4f   ',Best_FCAPSO_match_uniweight(repeat_times,i));
% %     end
% %     fprintf(fid,'\n');
% % end
% % fprintf(fid,'\n');
% % ave_Best_FCAPSO_match_uniweight = mean(Best_FCAPSO_match_uniweight);
% % std_Best_FCAPSO_match_uniweight = std(Best_FCAPSO_match_uniweight);
% % for i = 1:max_i
% %     fprintf(fid,'%12.4f   ',ave_Best_FCAPSO_match_uniweight(i));
% % end
% % fprintf(fid,'\n');
% % for i = 1:max_i
% %     fprintf(fid,'%12.4f   ',std_Best_FCAPSO_match_uniweight(i));
% % end
% % fclose(fid);
% % 
% fid = fopen('n_dFCAPSO_match_30.txt','wt');
% for repeat_times = 1:max_repeat    
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',Best_FCAPSO_match_divweight(repeat_times,i));
%     end
%     fprintf(fid,'\n');
% end
% fprintf(fid,'\n');
% ave_Best_FCAPSO_match_divweight = mean(Best_FCAPSO_match_divweight);
% std_Best_FCAPSO_match_divweight = std(Best_FCAPSO_match_divweight);
% for i = 1:max_i
%     fprintf(fid,'%12.4f   ',ave_Best_FCAPSO_match_divweight(i));
% end
% fprintf(fid,'\n');
% for i = 1:max_i
%     fprintf(fid,'%12.4f   ',std_Best_FCAPSO_match_divweight(i));
% end
% fclose(fid);   

% %%%FCAPSO 
% for i = 1:max_i
%     sigma_noise = 2*(i-1)/10;
%     for repeat_times = 1 : max_repeat
%         %t = cputime;
%         Best_FCAPSO_uniweight(repeat_times,i) = PSOed100(dim_obj,max_iteration,1,'FCAPSO',2);   %2 average
%         %t_FCAPSOu(i) = cputime - t;
%         %t = cputime;
%         %Best_FCAPSO_divweight(repeat_times,i) = PSOed100(dim_obj,max_iteration,2,'FCAPSO',2); 
%         %t_FCAPSOd(i) = cputime - t;        
%     end
% end
% 
% fid = fopen('n_uFCAPSO_30.txt','wt');
% for repeat_times = 1:max_repeat    
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',Best_FCAPSO_uniweight(repeat_times,i));
%     end
%     fprintf(fid,'\n');
% end
% fprintf(fid,'\n');
% ave_Best_FCAPSO_uniweight = mean(Best_FCAPSO_uniweight);
% std_Best_FCAPSO_uniweight = std(Best_FCAPSO_uniweight);
% for i = 1:max_i
%     fprintf(fid,'%12.4f   ',ave_Best_FCAPSO_uniweight(i));
% end
% fprintf(fid,'\n');
% for i = 1:max_i
%     fprintf(fid,'%12.4f   ',std_Best_FCAPSO_uniweight(i));
% end
% fclose(fid);  
% 
% fid = fopen('n_dFCAPSO_30.txt','wt');
% for repeat_times = 1:max_repeat    
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',Best_FCAPSO_divweight(repeat_times,i));
%     end
%     fprintf(fid,'\n');
% end
% fprintf(fid,'\n');
% ave_Best_FCAPSO_divweight = mean(Best_FCAPSO_divweight);
% std_Best_FCAPSO_divweight = std(Best_FCAPSO_divweight);
% for i = 1:max_i
%     fprintf(fid,'%12.4f   ',ave_Best_FCAPSO_divweight(i));
% end
% fprintf(fid,'\n');
% for i = 1:max_i
%     fprintf(fid,'%12.4f   ',std_Best_FCAPSO_divweight(i));
% end
% fclose(fid); 
%%FPSO
for i = 1:max_i
    sigma_noise = 2*(i-1)/10;
    for repeat_times = 1 : max_repeat
        % 每一列表示一个尺度
        % 每一行表示一次重复
        %t = cputime;
        Best_FPSO_uniweight(repeat_times,i) = PSOed100(dim_obj,max_iteration,1,'FPSO',1); 
        %t_FPSOu(i) = cputime- t;
        %t = cputime;
        %Best_FPSO_divweight(repeat_times,i) = PSOed100(dim_obj,max_iteration,2,'FPSO',1); 
        %t_FPSOd(i) = cputime- t;        
    end
end
fid = fopen('n_uFPSO_30.txt','wt');
for repeat_times = 1:max_repeat    
    for i = 1:max_i
        fprintf(fid,'%12.4f   ',Best_FPSO_uniweight(repeat_times,i));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
ave_Best_FPSO_uniweight = mean(Best_FPSO_uniweight);
std_Best_FPSO_uniweight = std(Best_FPSO_uniweight);
for i = 1:max_i
    fprintf(fid,'%12.4f   ',ave_Best_FPSO_uniweight(i));
end
fprintf(fid,'\n');
for i = 1:max_i
    fprintf(fid,'%12.4f   ',std_Best_FPSO_uniweight(i));
end
fclose(fid);

% fid = fopen('n_dFPSO_30.txt','wt');
% for repeat_times = 1:max_repeat   
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',Best_FPSO_divweight(repeat_times,i));
%     end
%     fprintf(fid,'\n');
% end
% fprintf(fid,'\n');
% ave_Best_FPSO_divweight = mean(Best_FPSO_divweight);
% std_Best_FPSO_divweight = std(Best_FPSO_divweight);
% for i = 1:max_i
%     fprintf(fid,'%12.4f   ',ave_Best_FPSO_divweight(i));
% end
% fprintf(fid,'\n');
% for i = 1:max_i
%     fprintf(fid,'%12.4f   ',std_Best_FPSO_divweight(i));
% end
% fclose(fid);

%%LPSO
for i = 1:max_i
    sigma_noise = 2*(i-1)/10;
    for repeat_times = 1 : max_repeat
        % 每一列表示一个尺度
        % 每一行表示一次重复
        %t = cputime;
        Best_LPSO_uniweight(repeat_times,i) = PSOed100(dim_obj,max_iteration,1,'LPSO',1);    %SPSO,FPSO,(FCAPSO,1)
        %t_LPSOu(i) = cputime- t;
        %t = cputime;
        %Best_LPSO_divweight(repeat_times,i) = PSOed100(dim_obj,max_iteration,2,'LPSO',1); 
        %t_LPSOd(i) = cputime- t;            
    end
end

fid = fopen('n_uLPSO_30.txt','wt');
for repeat_times = 1:max_repeat    
    for i = 1:max_i
        fprintf(fid,'%12.4f   ',Best_LPSO_uniweight(repeat_times,i));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
ave_Best_LPSO_uniweight = mean(Best_LPSO_uniweight);
std_Best_LPSO_uniweight = std(Best_LPSO_uniweight);
for i = 1:max_i
    fprintf(fid,'%12.4f   ',ave_Best_LPSO_uniweight(i));
end
fprintf(fid,'\n');
for i = 1:max_i
    fprintf(fid,'%12.4f   ',std_Best_LPSO_uniweight(i));
end
fclose(fid);
% fid = fopen('n_dLPSO_30.txt','wt');
% for repeat_times = 1:max_repeat    
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',Best_LPSO_divweight(repeat_times,i));
%     end
%     fprintf(fid,'\n');
% end
% fprintf(fid,'\n');
% ave_Best_LPSO_divweight = mean(Best_LPSO_divweight);
% std_Best_LPSO_divweight = std(Best_LPSO_divweight);
% for i = 1:max_i
%     fprintf(fid,'%12.4f   ',ave_Best_LPSO_divweight(i));
% end
% fprintf(fid,'\n');
% for i = 1:max_i
%     fprintf(fid,'%12.4f   ',std_Best_LPSO_divweight(i));
% end
% fclose(fid);

 
%%写文件

% fid = fopen('ProcessTime.txt','wt');
% for repeat_times = 1:max_repeat    
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',t_LPSOu(i));
%     end
%     fprintf(fid,'\n');
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',t_LPSOd(i));
%     end
%     fprintf(fid,'\n');
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',t_FPSOu(i));
%     end
%     fprintf(fid,'\n');
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',t_FPSOd(i));
%     end
%     fprintf(fid,'\n');
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',t_FCAPSOmu(i));
%     end
%     fprintf(fid,'\n');
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',t_FCAPSOmd(i));
%     end
%     fprintf(fid,'\n');
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',t_FCAPSOu(i));
%     end
%     fprintf(fid,'\n');
%     for i = 1:max_i
%         fprintf(fid,'%12.4f   ',t_FCAPSOd(i));
%     end
%     fprintf(fid,'\n');
% end
% fclose(fid);

% figure('Name','Time Performance')
% x_ax = 20:60:200;
% plot(x_ax,t_LPSOu,'b-+',x_ax,t_LPSOd,'b-x',...
%     x_ax,t_FPSOu,'g-o',x_ax,t_FPSOd,'g-*',...
%     x_ax,t_FCAPSOmu,'r-s',x_ax,t_FCAPSOmd,'r-d',...
%     x_ax,t_FCAPSOu,'m-p',x_ax,t_FCAPSOd,'m-h');
% legend('LPSOu','LPSOd','FPSOu','FPSOd','FCAPSOum','FCAPSOdm','FCAPSOu','FCAPSOd');

% figure('Name','Scalability Performance')
% x_ax = [2,50,100];
% plot(x_ax,ave_Best_LPSO_uniweight,'b-+',x_ax,ave_Best_LPSO_divweight,'b-x',...
%     x_ax,ave_Best_FPSO_uniweight,'g-o',x_ax,ave_Best_FPSO_divweight,'g-*',...
%     x_ax,ave_Best_FCAPSO_match_uniweight,'r-s',x_ax,ave_Best_FCAPSO_match_divweight,'r-d',...
%     x_ax,ave_Best_FCAPSO_uniweight,'m-p',x_ax,ave_Best_FCAPSO_divweight,'m-h');
% legend('LPSOu','LPSOd','FPSOu','FPSOd','FCAPSOum','FCAPSOdm','FCAPSOu','FCAPSOd');
% figure('Name','Noisy Performance')
% x_ax = 0:0.2:1.0;
% plot(x_ax,ave_Best_LPSO_uniweight,'b-+',x_ax,ave_Best_FPSO_uniweight,'g-o',...
%     x_ax,ave_Best_FCAPSO_match_divweight,'r-d');
% legend('LPSOu','FPSOu','FCAPSOdm');
figure('Name','Noisy Performance')
x_ax = 0:0.2:1.0;
plot(x_ax,ave_Best_LPSO_uniweight,'b-+',x_ax,ave_Best_FPSO_uniweight,'g-o');
legend('LPSOu','FPSOu');
%恒定多样权重
% i = 1;
% for dim_obj =   2:2 %10 : 10 : 100
%     Best_PSO_divweight(i) = PSOed100(dim_obj,1000,2,'FPSO');    
%     i = i + 1;
% end
% figure(2)
% plot(10:10:100, -Best_PSO_divweight,'b-+')
% legend('divweight');
% hold on


% 方差分析ANOVA
% 单因素方差分析 anoval
% i_group = 1;
% for dim_obj =  10 : 10 : 100
%     for repet = 1 : 12
%         Best_PSO(repet,i_group) = PSOed100(dim_obj);           
%     end
%     i_group = i_group +1;
% end
% p = anova1 (Best_PSO);
% 
% Best_PSO_uniweight = Best_PSO;
% p_uniweight = p;


% for repet = 1 : 1
%     Best_PSO(repet) = PSOed100(10,1000);   
%     %hold off
% end
% %figure(2)
% plot(Best_PSO);
% title(['Best PSO : Mean: ', num2str(mean(Best_PSO)), 'Std: ',...
%     num2str(std(Best_PSO))]);
