hold all
for dt = [.4 .2 .1 .05 .01]
   [t,g] = solve_for_g(2,dt);
   plot(t,g,'-*')
end
t = linspace(0,2);
g = .5 * exp(-t.^2);
plot(t,g,'-xk')
xlabel('t')
ylabel('y')
legend('dt = .4', 'dt = .2', 'dt = .1', 'dt = .05', 'dt = .01','Analytic soln')




function [tvec,g] = solve_for_g(T,dt)
%T is the total time to solve over, ie, (0,T)
%dt is step size
%initial condition is hard-coded in
g0 = .5;

N = ceil(T/dt);
tvec = 0:dt:N*dt;
g = zeros(1,N+1); %value of g at index i corresponds to g_{i-1}, ie,
% g = [g_0 g_1 g_2 ... g_N]
g(1) = g0;
for i = 0:N-1
    %matlab starts indexing at 1, so have to shift by 1 here
    g(i+2) = g(i+1) * (1-2*i*dt^2);  
end


end