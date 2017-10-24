hold all
for dt = [.4 .2 .1 .05 .01, .001]
   [t,g] = solve_for_g(4,dt);
   plot(t,g,'-o')
end
xlabel('x')
ylabel('y')
legend('dx = .4', 'dx = .2', 'dx = .1', 'dx = .05', 'dx = .01', 'dx = .001' )




function [tvec,g] = solve_for_g(T,dt)
%T is the total time to solve over, ie, (0,T)
%dt is step size
%initial condition is hard-coded in
g0 = 2;

N = ceil(T/dt);
tvec = 0:dt:N*dt;
g = zeros(1,N+1); %value of g at index i corresponds to g_{i-1}, ie,
% g = [g_0 g_1 g_2 ... g_N]
g(1) = g0;
for i = 1:N
    %matlab starts indexing at 1, so have to shift by 1 here
    num = g(i) + 2*dt*tvec(i+1) * cos(20*tvec(i+1));
    denom = 1 + 4*dt;
    g(i+1) = num/denom;
end


end