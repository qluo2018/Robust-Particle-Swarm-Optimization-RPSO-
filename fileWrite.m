function fileWrite(filename,data, max_repeat, max_i)
%%%Ð´ÎÄ¼þ
fid = fopen(filename,'wt');
for repeat_times = 1:max_repeat   
    for i = 1:max_i
        fprintf(fid,'%12.4f   ',data(repeat_times,i));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
ave_Best = mean(data);
std_Best = std(data);
for i = 1:max_i
    fprintf(fid,'%12.4f   ',ave_Best(i));
end
fprintf(fid,'\n');
for i = 1:max_i
    fprintf(fid,'%12.4f   ',std_Best(i));
end
fclose(fid);