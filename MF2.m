%%membership functions 
%%input: para����������ģ�����ϵĲ�����3*2�ľ���
%%output: x��3��������ֵ���ͣ��У���
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


