%%��������Ȩ�صĶ��ַ���
%%LPSO: Ȩ�ش�max��min�𽥼�С

function W_Next = uniUpdateWeight( iter, itmax, size_par, dim_obj)

weight_max = 0.9;  
weight_min = 0.4;
w_nex = weight_max - (( weight_max - weight_min ) / itmax ) * iter;

W_Next = w_nex * ones(size_par,dim_obj);
