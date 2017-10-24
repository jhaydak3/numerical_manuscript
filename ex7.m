close all
clc
clear

y0 = 1;
yN = 0;

[x1,y1] = mysolver(.1,y0,yN);
[x2,y2] = mysolver(.01,y0,yN);
[x3,y3] = mysolver(.001,y0,yN);
[x4,y4] = mysolver(.0005,y0,yN);

hold all
plot(x1,y1,'*-','MarkerSize',3)
plot(x2,y2,'*-','MarkerSize',3)
plot(x3,y3,'*-','MarkerSize',3)
plot(x4,y4,'*-','MarkerSize',3)
xlabel('x')
ylabel('y')
legend('dx = .1','dx = .01','dx = .001','dx = .0001')


function [xn,yn] = mysolver(dx,y0,yN)
xn = 1:dx:4;
yn = zeros(size(xn));
N = length(xn);
A = zeros(N,N);
b = zeros(N,1);
for i=2:N-1
    A(i-1,i-1) = f2(xn(i-1),dx);
    A(i-1,i) = -f1(xn(i-1),dx);
    A(i-1,i+1) = 1;
end
A(N-1,1) = 1;
b(N-1) = y0;
A(N,N) = 1;
b(N) = yN;
yn = A\b;

end

function out = f1(x,dx)
num = -x^2 + 1 + 2/(1+x^2)/(dx^2);
denom = 1/(1+x^2)/(dx^2) - sin(x)/2/dx;
out = num/denom;
end

function out  = f2(x,dx)
num = 1/(1+x^2)/(dx^2) + sin(x)/2/dx;
denom = 1/(1+x^2)/(dx^2) - sin(x)/2/dx;
out = num / denom;
end