global beta lambda Pr sigma n Sc E delta

beta =0.5;
lambda=0.1; 
Pr = 7;
sigma = 5;
n = 0.5;
Sc =1;
E =5;
delta = 1;

etaend = 8;

u(1) = -1.13;
u(2) = -0.23;
u(3) =-1.83;
u(4) = -0.55;

x0 = [0; 1; u(1); 0 ; u(2); 1; u(3); 1; u(4)];

figure(1)
step = 0.1;
options = odeset('AbsTol',1e-5,'Maxorder', 4);
[t,x] = ode45(@vdp, [0 :step: etaend], x0, options);
plot(t,x(:,2));

axis([0 etaend, -0.4, 1 ])
xlabel('\eta')
ylabel('f''(\eta)')
title('Figure 3')
legend('\lambda=0.1, \sigma=5, Pr=7, n=0.5 ,\beta=0.5, Sc=1, \delta=1')

function dxdeta = vdp(eta,x)
global beta lambda Pr sigma n Sc E delta
dxdeta = [x(2);
    x(3);
    (x(2)^2-x(1)*x(3)-2*lambda*(x(4)-beta*x(1)*x(5))-2*beta*x(1)*x(2)*x(3))/(1-beta*x(1)^2/2);...
    x(5);
    (x(2)*x(4)-x(1)*x(5)+2*lambda*(x(2)+beta*(x(2)^2-x(1)*x(3)+x(4)^2))-2*beta*x(1)*x(2)*x(5))*(1-beta*x(1)^2/2);...
    x(7);
    -Pr*x(1)*x(7);
    x(9);
    Sc*(sigma*(1+delta*x(6))^n*exp(-E/(1+delta*x(6)))*x(8)-x(1)*x(9))];
end