%�����Ƿ��������Բ���
%����Ҫע�⣬���ܱ�֤����̬�ֲ�
%����Wilcoxon�Ⱥͼ���

for i = 1:4
    LPSO = Best_div_LPSO(:,i);
    FPSO = Best_div_FPSO(:,i);
    RPSO = Best_roulette_RPSO(:,i);
    [pL(i), hL(i)] = ranksum(LPSO,RPSO,0.05);
    [pF(i), hF(i)] = ranksum(FPSO,RPSO,0.05);
end

hL
hF


%�����㷨Ч������Ժ�����ά���仯
%���е����ط������
%kruskalwallis����

group = [ones(1,30),2*ones(1,30),3*ones(1,30),4*ones(1,30)];
XL = [ Best_div_LPSO(:,1)', Best_div_LPSO(:,2)', Best_div_LPSO(:,3)', Best_div_LPSO(:,4)'];
XF = [ Best_div_FPSO(:,1)', Best_div_FPSO(:,2)', Best_div_FPSO(:,3)', Best_div_FPSO(:,4)'];
XR = [ Best_roulette_RPSO(:,1)', Best_roulette_RPSO(:,2)', Best_roulette_RPSO(:,3)', Best_roulette_RPSO(:,4)'];



% kruskalwallis
[pL,tableL,statsL] = anova1(XL,group)
[pF,tableF,statsF] = anova1(XF,group)
[pR,tableR,statsR] = anova1(XR,group)

% 
% subplot(3,1,1)
% boxplot(XL,group)
% title('LPSO for Rosenbrock Function')
% 
% subplot(3,1,2)
% boxplot(XF,group)
% title('FPSO for Rosenbrock Function')
% 
% subplot(3,1,3)
% boxplot(XR,group)
% title('RPSO for Rosenbrock Function')

