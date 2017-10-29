dx = 1e-3;
x0 = 1; xN = 5;
F = [];
wrange = -20:20;
for w = wrange
    IC = [0; w];
    [x,y] = myRK4(x0,IC,xN,dx);
    F = [F y(end)];
end
plot(wrange,F-1,'*-','MarkerSize',3)
xlabel('w')
ylabel('F(w,5) - 1')



function [x,y] = myRK4(x0,IC,xN,dx)
x = x0:dx:xN;
numx = length(x);
v = zeros(2,numx);
v(1,1) = IC(1);
v(2,1) = IC(2);

for i = 1:numx-1
    vn = v(:,i);
    h1 = dx*A(x(i))*vn;
    h2 = dx*A(x(i)+.5*dx)*(vn+.5*h1);
    h3 = dx*A(x(i)+.5*dx)*(vn+.5*h2);
    h4 = dx*A(x(i+1))*(vn+h3);
    v(:,i+1) = vn + (1/6)*(h1+2*h2+2*h3+h4);
end
y = v(1,:);
end

function out =  A(x)
out = [0 1; (5/x^2-1) -1];
end