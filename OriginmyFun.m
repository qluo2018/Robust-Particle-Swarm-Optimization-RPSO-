%目标函数非负
function z = OriginmyFun(x)
z = 0;
n = size(x,2);
%y = 100*(x(2)-x(1)^2)^2+(x(1)-1)^2;
% %%Schaffer function 
% for i =1:n
%     z = z + x(i)^2;
% end 
% z = 0.5 + ((sin(sqrt(z)))^2 - 0.5) / (1 + 0.001*z)^2;

% % %%Beale function just for two dimensions
% y1 = 1.5; y2 = 2.25; y3 = 2.625;
% z = (y1 - x(1)*(1-x(2)))^2 + (y2 - x(1)*(1 - x(2)^2))^2 + (y3 - x(1)*(1-x(2)^3)^2;

% %%Rastrgian function
% for i =1 : n
%     z = z + x(i) * x(i) - 10 * cos(2*pi*x(i)) + 10;
% end  
%%Rossenbrock function 
% for i = 1 : n-1
%     z = z + 100 * (x(i+1) - x(i) * x(i)) * (x(i+1) - x(i) * x(i));
%     z = z + ( x(i) - 1 ) * ( x(i) - 1 );
% end
%% Six-humpCamelback %-50,500
%  for i = 1: n/2
%      z = z + (4+2.1*x(i)*x(i)+0.33333*x(i)*x(i)*x(i)*x(i))*x(i)*x(i);
%      z = z + x(i)*x(i+n/2) + (-4 + 4 * x(i+n/2)*x(i+n/2)*x(i+n/2)*x(i+n/2));
%  end
% % %% Sphere %0,500
for i = 1:n
    z = z + x(i)*x(i);
end

%% 2n minima
% for i = 1 : n
%     z = z + x(i)*x(i)*x(i)*x(i);
%     z = z + 5*x(i);
%     z = z - 16*x(i)*x(i);
% end

%%% Griewank
% z1 = 0;
% z2 = 1;
% for i = 1:n
%     z1 = z1 + x(i)*x(i);
%     z2 = z2 * cos(x(i)/sqrt(i));
% end
% z = z1/4000 - z2 + 1;
