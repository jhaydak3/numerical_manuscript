dx = 1e-2;
x0 = 0; xN = 10;
F = [];
wrange = -50:1:50;
beta = .5; alpha = 2;
for w = wrange
    IC = [alpha; w];
    [x,y] = myRK4(x0,IC,xN,dx);
    F = [F y(end)];
end
g = @(xq) interp1(wrange,F-beta,xq);
hold on
plot(wrange,F-beta,'*-','MarkerSize',3)
xlabel('w')
ylabel('F(w,10) - .5')


function [x,y] = myRK4(x0,IC,xN,dx)
x = x0:dx:xN;
numx = length(x);
v = zeros(2,numx);
v(1,1) = IC(1);
v(2,1) = IC(2);

for i = 1:numx-1
    vn = v(:,i);
    h1 = dx*dV(x(i),vn);
    h2 = dx*dV(x(i)+.5*dx,vn+.5*h1);
    h3 = dx*dV(x(i)+.5*dx,vn+.5*h2);
    h4 = dx*dV(x(i+1),vn+h3);
    v(:,i+1) = vn + (1/6)*(h1+2*h2+2*h3+h4);
end
y = v(1,:);
end

function out =  dV(x,v)
out = zeros(2,1);
out(1) = v(2);
out(2) = -x*v(1) -v(1)^2*v(2);
end