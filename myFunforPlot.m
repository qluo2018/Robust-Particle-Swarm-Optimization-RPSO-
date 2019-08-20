%目标函数图形显示
function myFunPlot(x,y)
% x = [-5:0.1:5];
% y = [-5:0.1:5];
% 
% for i = 1:101
%     for j = 1:101
%         z(i,j) = myFun([x(i),y(j)]);
%     end
% end
% 
% surfc(x,y,z)
% colormap hsv

%%Rosenbrock
%f = ['100 * ( x - y^2 )^2 + (x - 1)^2 '];

%f = ['x^2+y^2']; %Sphere%
%%Rastrigrin
  f = ['x^2 - 10 * cos(2*pi*x)'...
     '+y^2 - 10 * cos(2*pi*y)'...
     '+20'];
%f = ['(4+2.1*x^2+1/3*x^4)*x^2+x*y+(-4+4*y^2)*y']; %six-humpCamlback
%f = ['1/4000*(x^2+y^2) - cos(x)*cos(y/sqrt(2)) + 1'];

ezsurfc(f,[-5,5],[-5,5])   
%ezsurf(f,[-5,5],[-5,5])
%ezcontourf(f,[-5,5],[-5,5])
%ezcontour('x^2 - 10 * cos(2*pi*x)+ y^2 - 10 * cos(2*pi*y) + 20',[-5,5],[-5,5])
colormap hsv
% f = ['3*(1-x)^2*exp(-(x^2)-(y+1)^2)',... 
%      '- 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2)',... 
%      '- 1/3*exp(-(x+1)^2 - y^2)'];
%  ezcontourf(f,[-3,3],49)