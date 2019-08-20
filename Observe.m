%%观察粒子的移动与最优路径的差距
clear
clc
function_name = 'Rastrigian';
clear ObjPara;
ObjPara = struct('NCBP',[0, 100], 'Dim', 30, 'Name', function_name);
clear AlgPara;   
AlgPara = struct('MaxIte',30, 'WeightType', 1, ...
            'MatchType', 1, 'IniRange',[-5, +5],...
            'Algorithm', 'LPSO', 'Size', 20);
        
%  Find_uni_LPSO = PSOed200(ObjPara,AlgPara); 
% 
%  clear AlgPara.Algorithm;
%  AlgPara.Algorithm = 'FPSO';
%  Find_uni_FPSO = PSOed200(ObjPara,AlgPara); 
%  
%  
clear AlgPara.Algorithm;
AlgPara.Algorithm = 'RPSO';
AlgPara.WeightType = 2; 
Find_fixed_RPSO = PSOed200(ObjPara,AlgPara);  