%%��������Ȩ�صĶ��ַ���
%%LPSO: Ȩ�شӳ�ʼȨ�ص�min�𽥼�С

function W_Next = divUpdateWeight( iter, itmax, ini_W)

weight_min = 0.4;
W_Next = ini_W - (( ini_W - weight_min ) / itmax ) * iter;
