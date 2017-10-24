hold all
x0 = .75; y0=0; xn = 10;
[x1,y1] = myRK2(.25,x0,y0,xn);
[x2,y2] = myRK4(.25,x0,y0,xn);
[x3,y3] = my_jon(.25,x0,y0,xn);
[x4,y4] = myRK2(.05,x0,y0,xn);
[x5,y5] = myRK4(.05,x0,y0,xn);
[x6,y6] = my_jon(.05,x0,y0,xn);
[x_c,y_c] = myRK4(.001,x0,y0,xn);
plot(x1,y1,'r-o',x2,y2,'c-o',x3,y3,'b-o')
plot(x4,y4,'r-o','MarkerFaceColor','r','MarkerSize',3)
plot(x5,y5,'c-o', 'MarkerFaceColor','c','MarkerSize',3)
plot(x6,y6,'b-o', 'MarkerFaceColor','b','MarkerSize',3)
plot(x_c,y_c,'k-x','MarkerSize',2)
xlabel('x')
ylabel('y')
legend('dx = .25, RK2', 'dx = .25, RK4', 'dx = .25, jon', 'dx = .05, RK2', ...
    'dx = .05, RK4', 'dx = .05, jon','dx = .001, RK4')

function [x,y] = myRK2(dx,x0,y0,xn)
%dx - timestep, x0 - initial point, y0- initial point, xn - final x to be
%integrated over
x = x0:dx:xn;
y = zeros(size(x));
y(1) = y0;
numX = length(x);
for i = 2:numX
   h1 = dx * dydx(x(i-1),y(i-1));
   h2 = dx * dydx(x(i-1)+.5*dx,y(i-1)+.5*h1);
   y(i) = y(i-1) + h2;
end
end

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

function [x,y] = my_jon(dx,x0,y0,xn)
%dx - timestep, x0 - initial point, y0- initial point, xn - final x to be
%integrated over
x = x0:dx:xn;
y = zeros(size(x));
y(1) = y0;
numX = length(x);
for i = 2:numX
   h1 = dx * dydx(x(i-1),y(i-1));
   h2 = dx * dydx(x(i-1)+.25*dx,y(i-1)+.25*h1);
   h3 = dx * dydx(x(i-1)+.5*dx,y(i-1)+.5*h2);
   h4 = dx * dydx(x(i-1)+.75*dx,y(i-1)+.75*h3);
   h5 = dx * dydx(x(i),y(i-1)+h4);
   y(i) = y(i-1) + (1/9) * (h1 + 2*h2 + 3*h3 + 2*h4 + h5);
end
end

function out = dydx(x,y)
out = sin(x)/x^4 + 2/x^4 - y/x - 2*y/x^3;
end