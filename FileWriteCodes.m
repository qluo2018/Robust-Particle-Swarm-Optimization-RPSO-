%Ð´ÎÄ¼þµÄÓï¾ä
fid = fopen('r2_uLPSO_30.txt','wt');
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

fid = fopen('r2_dLPSO_30.txt','wt');
for repeat_times = 1:max_repeat    
    for i = 1:max_i
        fprintf(fid,'%12.4f   ',Best_LPSO_divweight(repeat_times,i));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
ave_Best_LPSO_divweight = mean(Best_LPSO_divweight);
std_Best_LPSO_divweight = std(Best_LPSO_divweight);
for i = 1:max_i
    fprintf(fid,'%12.4f   ',ave_Best_LPSO_divweight(i));
end
fprintf(fid,'\n');
for i = 1:max_i
    fprintf(fid,'%12.4f   ',std_Best_LPSO_divweight(i));
end
fclose(fid);

fid = fopen('r2_uFPSO_30.txt','wt');
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

fid = fopen('r2_dFPSO_30.txt','wt');
for repeat_times = 1:max_repeat   
    for i = 1:max_i
        fprintf(fid,'%12.4f   ',Best_FPSO_divweight(repeat_times,i));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
ave_Best_FPSO_divweight = mean(Best_FPSO_divweight);
std_Best_FPSO_divweight = std(Best_FPSO_divweight);
for i = 1:max_i
    fprintf(fid,'%12.4f   ',ave_Best_FPSO_divweight(i));
end
fprintf(fid,'\n');
for i = 1:max_i
    fprintf(fid,'%12.4f   ',std_Best_FPSO_divweight(i));
end
fclose(fid);

fid = fopen('r_uFCAPSO_match_30.txt','wt');
for repeat_times = 1:max_repeat   
    for i = 1:max_i
        fprintf(fid,'%12.4f   ',Best_FCAPSO_match_uniweight(repeat_times,i));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
ave_Best_FCAPSO_match_uniweight = mean(Best_FCAPSO_match_uniweight);
std_Best_FCAPSO_match_uniweight = std(Best_FCAPSO_match_uniweight);
for i = 1:max_i
    fprintf(fid,'%12.4f   ',ave_Best_FCAPSO_match_uniweight(i));
end
fprintf(fid,'\n');
for i = 1:max_i
    fprintf(fid,'%12.4f   ',std_Best_FCAPSO_match_uniweight(i));
end
fclose(fid);

fid = fopen('r_dFCAPSO_match_30.txt','wt');
for repeat_times = 1:max_repeat    
    for i = 1:max_i
        fprintf(fid,'%12.4f   ',Best_FCAPSO_match_divweight(repeat_times,i));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
ave_Best_FCAPSO_match_divweight = mean(Best_FCAPSO_match_divweight);
std_Best_FCAPSO_match_divweight = std(Best_FCAPSO_match_divweight);
for i = 1:max_i
    fprintf(fid,'%12.4f   ',ave_Best_FCAPSO_match_divweight(i));
end
fprintf(fid,'\n');
for i = 1:max_i
    fprintf(fid,'%12.4f   ',std_Best_FCAPSO_match_divweight(i));
end
fclose(fid);      

fid = fopen('r_uFCAPSO_30.txt','wt');
for repeat_times = 1:max_repeat    
    for i = 1:max_i
        fprintf(fid,'%12.4f   ',Best_FCAPSO_uniweight(repeat_times,i));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
ave_Best_FCAPSO_uniweight = mean(Best_FCAPSO_uniweight);
std_Best_FCAPSO_uniweight = std(Best_FCAPSO_uniweight);
for i = 1:max_i
    fprintf(fid,'%12.4f   ',ave_Best_FCAPSO_uniweight(i));
end
fprintf(fid,'\n');
for i = 1:max_i
    fprintf(fid,'%12.4f   ',std_Best_FCAPSO_uniweight(i));
end
fclose(fid);  

fid = fopen('r_dFCAPSO_30.txt','wt');
for repeat_times = 1:max_repeat    
    for i = 1:max_i
        fprintf(fid,'%12.4f   ',Best_FCAPSO_divweight(repeat_times,i));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
ave_Best_FCAPSO_divweight = mean(Best_FCAPSO_divweight);
std_Best_FCAPSO_divweight = std(Best_FCAPSO_divweight);
for i = 1:max_i
    fprintf(fid,'%12.4f   ',ave_Best_FCAPSO_divweight(i));
end
fprintf(fid,'\n');
for i = 1:max_i
    fprintf(fid,'%12.4f   ',std_Best_FCAPSO_divweight(i));
end
fclose(fid);  

% Result = [Best_LPSO_uniweight;Best_LPSO_divweight;Best_FPSO_uniweight;Best_FPSO_divweight;...
%  Best_FCAPSO_match_uniweight;Best_FCAPSO_match_divweight;Best_FCAPSO_uniweight;Best_FCAPSO_divweight];
