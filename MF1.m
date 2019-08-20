%%membership functions 
%%input: para包含了三个模糊集合的参数，3*2的矩阵
%%output: x的3个隶属度值：低，中，高
function mu_Low = MF1 ( x, para )


%Low
mu_Low = 0;
if x < para(1,1)
    mu_Low = 1;
elseif x > para(1,2)
    mu_Low = 0;
else
    mu_Low = ( para(1,2) - x ) / ( para(1,2) - para(1,1) );
end

%Medium
mu_Med = 0;
if x < para(2,1)
    mu_Med = 0;
elseif x > para(2,2)
    mu_Med = 0;
elseif ( x > para(2,1) && x < ( para(2,1) + para(2,2) ) * 0.5 )
    mu_Med = 2 * ( x - para(2,1) ) / ( para(2,2) - para(2,1) );   
else
    mu_Med = 2 * ( para(2,2) - x ) / ( para(2,2) - para(2,1) );   
end

%High
mu_Hig = 0;
if x < para(3,1)
    mu_Hig = 0;
elseif x > para(3,2)
    mu_Hig = 1;
else
    mu_Hig = (  x - para(3,1) ) / ( para(3,2) - para(3,1) );
end
