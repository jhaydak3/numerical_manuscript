hold all
[x1,y1] = solve_for_y(10,.4);
[x2,y2] = solve_for_y(10,.2);
[x3,y3] = solve_for_y(10,.1);
[x4,y4] = solve_for_y(10,.05);
[x5,y5] = solve_for_y(10,.01);
[x6,y6] = solve_for_y(10,.001);
plot(x1,y1,'-o',x2,y2,'-o',x3,y3,'-o',x4,y4,'-',x5,y5,'-',x6,y6,'-')
xlabel('x')
ylabel('y')
legend('dx = .4', 'dx = .2', 'dx = .1', 'dx = .05', 'dx = .01', 'dx = .001' )
axis([0 10 .8 2.6])
function [xvec,y] = solve_for_y(T,dx)
%T is the total time to solve over, ie, (0,T)
%dx is step size
%initial condition is hard-coded in
y0 = 2;

N = ceil(T/dx);
xvec = 0:dx:N*dx;
y = zeros(1,N+1); %value of y at index i corresponds to y_{i-1}, ie,
% y = [y_0 y_1 y_2 ... y_N]
y(1) = y0;
for i = 1:N
    %matlab starts indexing at 1, so have to shift by 1 here
    y_p = y(i) + dx * (2*xvec(i)*cos(25*xvec(i))^2 - y(i)^3);
    y(i+1) = y(i) + (dx/2) * (2*xvec(i+1)*cos(25*xvec(i+1))^2 - y_p^3 + 2*xvec(i) * cos(25*xvec(i))^2 - y(i)^3);
end
end