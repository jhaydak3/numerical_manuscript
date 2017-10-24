
[x1,y1] = myfun(1);
[x2,y2] = myfun(.1);
[x3,y3] = myfun(.01);
[x4,y4] = myfun(.001);

figure
hold all
plot(x1,y1,'*-',x2,y2,'^-',x3,y3,'o-',x4,y4,'x-')
legend('dx = 1', 'dx = .1', 'dx = .01', 'dx = .001')
xlabel('x')
ylabel('y')
axis([1 10 -inf inf])

function [x,y] =  myfun(dx)
x = 1:dx:(10+2*dx); %have to go a little over to be able to solve for y(10)
y = x; %get correct size of y

y(1) = .75; %y_0 = .75
y(2) = y(1) + dx; %approx y2 with euler method

num_x = length(x);
for i = 3:num_x
    num= 4*dx^2/(x(i-1)^2) - dx^2 + 2;
    denom = 1 + dx/(2*x(i-1));
    A = num/denom;
    num = dx/(2*x(i-1)) - 1;
    B = num/denom;
    y(i) = A*y(i-1)+ B*y(i-2);
end
end