%% ����GA �Ի����������Ŵ�����
% function GeneticOperations ( Population )
% 
% %%----�Ŵ�����������������----%%
% CrossoverProbability = 0.80;
% MutationProbability = 0.001;
% MaximumNumberOfGenerations =  10;
% 
% 
% [x fval] = ga(@rastriginsfcn, 2)


options = gaoptimset('Generations',300);
rand('state', 71); % These two commands are only included to
randn('state', 59); % make the results reproducible.
record=[];
for n=0:.05:1
 options = gaoptimset(options,'CrossoverFraction', n);
 [x fval]=ga(@rastriginsfcn, 10,[],[],[],[],[],[],[],options);
 record = [record; fval];
end
plot(0:.05:1, record);
xlabel('Crossover Fraction');
ylabel('fval')