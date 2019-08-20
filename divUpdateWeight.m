%%调整惯性权重的多种方法
%%LPSO: 权重从初始权重到min逐渐减小

function W_Next = divUpdateWeight( iter, itmax, ini_W)

weight_min = 0.4;
W_Next = ini_W - (( ini_W - weight_min ) / itmax ) * iter;
