close all
clc
clear
dx_init = .01;
dx_max = .5;
tol = 1e-6;
x0 =2; y0 = 2;
x_final = 10;
thisx = x0; thisy = y0; thisdx = .01;
k  = 5; %one iteration of RK4 is 5th order accurate
xvec = [x0];
yvec = [y0];
while thisx < x_final
    %calculate A(h)
    Ah = RK4_step(thisdx,thisx,thisy);
    Ah2 = RK4_step(thisdx/2,thisx,thisy);
    Ah2 = RK4_step(thisdx/2,thisx+thisdx/2,Ah2);
    B = (Ah - Ah2) / thisdx^k;
    B = B/(1 - (1/2^k));
    Err = B * (thisdx/2)^k;
    if abs(Err) <= tol %then accept
        thisy = ((2^k)*Ah2-Ah)/(2^k-1);
        thisx = thisx + thisdx;
        xvec = [xvec thisx];
        yvec = [yvec thisy];
    end
    thisdx = abs((tol*10/B))^(1/k);
    %thisdx = thisdx*.9;
end
hold all
plot(xvec,yvec,'*')
[x1,y1] = myRK4(.01,x0,y0,x_final);
plot(x1,y1,'-')
xlabel('x')
ylabel('y')
legend('Adaptive RK4','Uniform RK4, dx = .01')

function [x,y] = myRK4(dx,x0,y0,xn)
%dx - timestep, x0 - initial point, y0- initial point, xn - final x to be
%integrated over
x = x0:dx:xn;
y = zeros(size(x));
y(1) = y0;
numX = length(x);
for i = 2:numX
   h1 = dx * dydx(x(i-1),y(i-1));
   h2 = dx * dydx(x(i-1)+.5*dx,y(i-1)+.5*h1);
   h3 = dx * dydx(x(i-1)+.5*dx,y(i-1)+.5*h2);
   h4 = dx * dydx(x(i),y(i-1)+h3);
   y(i) = y(i-1) + (1/6) * (h1 + 2*h2 + 2*h3 + h4);
end
end


function [yn] = RK4_step(dx,x0,y0)
%dx - timestep, x0 - initial point, y0- initial point, xn - final x to be
%integrated over
h1 = dx * dydx(x0,y0);
h2 = dx * dydx(x0+.5*dx,y0+.5*h1);
h3 = dx * dydx(x0+.5*dx,y0+.5*h2);
h4 = dx * dydx(x0+dx,y0+h3);
yn = y0 + (1/6) * (h1 + 2*h2 + 2*h3 + h4);
end

function out = dydx(x,y)
out = (.01*x^2-2)*y^.5 + exp(-x^2)*y+x^2*sin(x)^2;
end