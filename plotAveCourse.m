%画平均迭代过程
function plotAveCourse(BestGA, AveGA, BestPSO, AvePSO, Name, Title1, Title2, i)

    GAFitness  =  mean( BestGA,2);
    meanGAFitness   =  mean( AveGA,2);
    GBest =  mean( BestPSO,2);
    meanLBest  =  mean( AvePSO,2);   
    Name = [Name, ' for ',num2str(20 + 60 * (i-1))];
    figure('Name', Name);
    x = 1: size(GBest,1);
    subplot(2,1,1)
    plot(x, GAFitness,'r-*',x,meanGAFitness,'b-+');
    title(Title1)
    subplot(2,1,2)
    plot(x,GBest,'r-*',x,meanLBest,'b-+')    
    title(Title2)
    
    hgsave(Name)