global beta lambda Pr sigma n Sc E delta

beta =0.5;
lambda=0.2; 
Pr = 7;
sigma = 1;
n = 0.5;
Sc =1;
E =5;
delta = 1;
etaend = 10;

u(1) = -1.13;
u(2) = -0.23;
u(3) =-1.84545;
u(4) = -0.55;

ue(1) = 0;
ue(2) = 0.5;
ue(3) = -0.51903;
ue(4) = -0.95064;

x0 = [0; 1; u(1); 0 ; u(2); 1; u(3); 1; u(4)];

step = 0.1;

options1 = odeset('AbsTol',1e-2,'Maxorder', 4);
options2 = odeset('AbsTol',1e-2,'Maxorder', 2);

[t,x1] = ode45(@vdp, [0 :step: etaend], x0, options1);
[t,x2] = ode113(@vdp, [0 :step: etaend], x0, options2);
[t,x3] = euler(@vdp, x0, etaend, step);

plot(t,x1(:,6),'r',t,x2(:,6),'g',t,x3(6,:),'b');
axis([0 etaend 0 1])
xlabel('\eta')
ylabel('\theta(\eta)')
title('Three methods')
legend('RK4','Adams','Euler')

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