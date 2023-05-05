function z=Fshoot(u)

x0 = [0; 1; u(1); 0 ; u(2); 1; u(3); 1; u(4)];
[t,x] = ode45(@odefun20, 0:0.001:5, x0);

z=[x(end,2);x(end,4); x(end,6);x(end,8)];
