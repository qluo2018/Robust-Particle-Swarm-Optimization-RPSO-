function y = ComputeObjFunctionForParticle(x, Name)

y = [0];
[num_par, dim_obj, it] = size(x);
for i = 1:num_par
    x_temp(1,:) = x(i,:,1);
    y(i,1) = myFun(x_temp, Name);    
end

