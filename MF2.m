%%membership functions 
%%input: para包含了三个模糊集合的参数，3*2的矩阵
%%output: x的3个隶属度值：低，中，高
function mu_Med = MF2 ( x, para )


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


