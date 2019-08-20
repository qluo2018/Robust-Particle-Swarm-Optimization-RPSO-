%目标函数非负
function z = myFun(x,Name)

z = 0;
n = size(x,2);

switch Name
    case  'Schaffer'
        % 有无数个极小值点，其中只有一个（0，0）为最小值点，最小值为0
        % 取值范围[-100,100]^2,最小值周围一圈值都是0.0097
        for i =1:n
            z = z + x(i)^2;
        end
        z = 0.5 + ((sin(sqrt(z)))^2 - 0.5) / (1 + 0.001*z)^2;
    case 'Rastrigian'
        % 最小值在（0，0）取到0
        for i =1 : n
            z = z + x(i) * x(i) - 10 * cos(2*pi*x(i)) + 10;
        end       
    case 'Rossenbrock' %DeJone function 2
        % 最小值在（0，0）取到0
        for i = 1 : n-1
            z = z + 100 * (x(i+1) - x(i) * x(i)) * (x(i+1) - x(i) * x(i));
            z = z + ( x(i) - 1 ) * ( x(i) - 1 );
        end
    case 'Six-humpCamelback' 
        %2维的时候有6个极小值点，其中(-0.0898,0.7126)和(0.0898,-0.7126)为最小值点，
        %最小值为-1.031628，自变量的范围[-5，5]^2
        % -50
        for i = 1: n/2
            z = z + (4+2.1*x(i)*x(i)+0.33333*x(i)*x(i)*x(i)*x(i))*x(i)*x(i);
            z = z + x(i)*x(i+n/2) + (-4 + 4 * x(i+n/2)*x(i+n/2)*x(i+n/2)*x(i+n/2));
        end
    case 'Sphere'  %DeJone function 1
        %最小值在0处取到0
        for i = 1:n
            z = z + x(i)*x(i);
        end
    case 'Weierstrass' 
        %%最小值在0处取到0
        a = 0.5; b = 3; kmax = 20;
        z = 0;
        for i = 1 : n
            s1 = 0;            
            for k = 1 : kmax
                s1 = s1 + a^k*cos(2*pi*b^k*(x(i)+0.5));                
            end
            z = z + s1;
        end
        s2 = 0;
        for k = 1 : kmax
            s2 = s2 + a^k*cos(2*pi*b^k*0.5);
        end
        z = z - n * s2; 
    case '2n-minima'
        % -15000~-20000不清楚
        for i = 1 : n
            z = z + x(i)*x(i)*x(i)*x(i);
            z = z + 5*x(i);
            z = z - 16*x(i)*x(i);
        end
    case 'Griewank'
        %最小值在0处取到，0
        z1 = 0;z2 = 1;
        for i = 1:n
            z1 = z1 + x(i)*x(i);
            z2 = z2 * cos(x(i)/sqrt(i));
        end
        z = z1/4000 - z2 + 1;
    case 'Step'  %%De Jone 3 最小值0
        z = 5 * n;
        for i = 1 : n
            z = z + floor(x(i) + 0.5) ^2;            
        end
    case 'Schwefels-Problem-2.22'  %% 最小值为0
        z = 0;
        s = 1;
        for i = 1: n
            z = z + abs(x(i));
            s = s * abs(x(i));
        end
        z = z + s;
    case 'Schwefels-Problem-1.2'  %%最小值0
        z = 0;
        for i = 1 : n
            s = 0;
            for j = 1 : i
                s = s + x(i);
            end
            z = z + s*s;
        end
    case 'Schwefels-Problem-2.21' %%最小值0
        z = max(abs(x));
    case 'Quartic'  %%最小值0
        z = 0;
        for i = 1:n
            z = z + i * x(i)^4;
        end
        z = z  + random('Normal',0,1);
    case 'Schwefels-Problem-2.26' %%30维的时候, f8（ x*） = f8 （ 420. 9687，?，420.9687）= - 12569. 5
        z = 0;
        for i = 1 :n
            z = z + x(i) * sin(sqrt(abs(x(i))));
        end
        z = -z;
    case 'Schwefel'
        %%最小值[420.96， 。。。， 420.96]取0
        z = 0;
        for i = 1 : n
            z = z + x(i) * sin(sqrt(abs(x(i))));
        end
        z  = 418.9829 * n - z;
    case 'Ackley' %%最小值为0
        z = 0;
        s = 0;
        s1 = 0;
        for i = 1 : n
            s = s + x(i) *x(i);
            s1 = s1 + cos(2*pi*x(i));
        end
        z = -20 * exp(-0.2*sqrt(s/n)) - exp(s1/n) + 20 + exp(1);
%     case 'Penalized' %%最小值在1取到0
%         z = 0;
%         s = 0;
%         for i = 1 : n-1
%             s = s + (x(i)-1)*(x(i)-1)*(1+sin(3*pi*x(i+1)));
%         end
%         z = 0.1 * (sin(3*pi*x(1)) + s + (x(n)-1)*(x(n)-1));
%         s = 0;
%         for i = 1 : n
%             s = s + 
%             
    case 'Shekel'
        z = 0;
        a = [-32, -16, 0, 16, 32, -32, -16, 0, 16, 32,-32, -16, 0, 16, 32,...
            -32, -16, 0, 16, 32,-32, -16, 0, 16, 32;
            -32, -32, -32,-32,-32,-16,-16,-16,-16,-16,...
            0,0,0,0,0,16,16,16,16,16, 32, 32, 32, 32, 32];
        s2 = 0;
        for j = 1:25 
            s = 0;
            for i = 1: n
                s = s + (x(i) - a(i,j))^6;
            end
            s2 = 1.0 / (j + s);
        end
        z = s2 + 1/500;
        z = 1 / z; 
end

% % SPHERE is the dream of every optimization algorithm. It is smooth, it is
% unimodal, it is symmetric and it does not have any of the problems we have
% discussed so far. The performance on SPHERE is a measure of the general 
% efficiency of the algorithm. 
% % ROSENBROCK is the nightmare. It has a very narrow ridge. The tip of the
% ridge is very sharp, and it runs around a parabola. Algorithms that are 
% not able to discover good directions underperform in this problem. 
% %STEP is the representative of the problem of flat surfaces. Flat surfaces 
% are obstacles for optimization algorithms, because they do not give any
% information as to which direction is favorable. Unless an algorithm has 
% variable step sizes, it can get stuck on one of the flat plateaus. 
% % QUARTIC is a simple unimodal function padded with noise. The gaussian 
% noise makes sure that the algorithm never gets the same value on the same
% point. Algorithms that do not do well on this test function will do poorly 
% on noisy data. 
% % FOXHOLES is an example of many (in this case 25) local optima. Many 
% standard optimization algorithms get stuck in the first peak they find.    
         
        