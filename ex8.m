close all
x0 = 0;
v0 = [0;6;-5.5];
xn = 12;
[x1,v1] = my_RK4(.01,x0,v0,xn);
[x2,v2] = my_RK4(.005,x0,v0,xn);
[x3,v3] = my_RK4(.001,x0,v0,xn);
hold all
plot(x1,v1,'-o','MarkerSize',3)
plot(x2,v2,'-^','MarkerSize',3)
plot(x3,v3,'-','MarkerSize',3)
xlabel('x')
ylabel('v')
legend('dx = .01','dx = .005','dx = .001')
function [xvec, v] = my_RK4(dx,x0,v0,xn)
%dx - timestep, x0 - initial point, v0- initial point, xn final time
xvec = x0:dx:xn;
v = [v0(1)];
for x0 = xvec(1:end-1)
    h1 = dx * dvdx(x0,v0);
    h2 = dx * dvdx(x0+.5*dx,v0+.5*h1);
    h3 = dx * dvdx(x0+.5*dx,v0+.5*h2);
    h4 = dx * dvdx(x0+dx,v0+h3);
    v0 = v0 + (1/6) * (h1 + 2*h2 + 2*h3 + h4); 
    v = [v v0(1)];
end

end

function out = dvdx(x,v)
out = zeros(3,1);
out(1) = v(2);
out(2) = v(3);
out(3) = -2*x^2 * v(3) - x*v(2);
end