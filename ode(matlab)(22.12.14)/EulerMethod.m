global beta lambda Pr sigma n Sc E delta
% beta = 0;
% Pr = 7;
% sigma = 5;
% n = 0.5;
% Sc = .25;
% E = 5;
% delta = 1;
% etaend = 4;
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


%% effect accoring ramda
lambda1 = 1;

figure(1)
step_v = [0.1;0.5;1];
for i=1:length(step_v)
 step=step_v(i);
 [t,x] = euler(@vdp, x0, etaend, step);
 plot(t,x(2,:));
 hold on
end
hold off
axis([0 etaend -0.2 1])
xlabel('\eta')
ylabel('f''(\eta)')
title('Euler')
legend('step =0.1','step =0.5','step =1')

figure(2)
for i=1:length(step_v)
 step=step_v(i);
 [t,x] = euler(@vdp, x0, etaend, step);
 plot(t,x(4,:));
 hold on
end
hold off
axis([0 etaend -1 1])
xlabel('\eta')
ylabel('g(\eta)')
title('Euler')
legend('step =0.1','step =0.5','step =1')

figure(3)
for i=1:length(step_v)
 step=step_v(i);
 [t,x] = euler(@vdp, x0, etaend, step);
 plot(t,x(6,:));
 hold on
end
hold off
axis([0 etaend 0 1])
xlabel('\eta')
ylabel('\theta(\eta)')
title('Euler')
legend('step =0.1','step =0.5','step =1')

figure(4)
for i=1:length(step_v)
 step=step_v(i);
 [t,x] = euler(@vdp, x0, etaend, step);
 plot(t,x(8,:));
 hold on
end
hold off
axis([0 etaend 0 1])
xlabel('\eta')
ylabel('\phi(eta)')
title('Euler')
legend('step =0.1','step =0.5','step =1')

function dxdeta = vdp(eta,x)
global beta lambda Pr sigma n Sc E delta
dxdeta = [x(2)
    x(3)
    (x(2)^2-x(1)*x(3)-2*lambda*(x(4)-beta*x(1)*x(5))-2*beta*x(1)*x(2)*x(3))/(1-beta*x(1)^2/2)
    x(5)
    (x(2)*x(4)-x(1)*x(5)+2*lambda*(x(2)+beta*(x(2)^2-x(1)*x(3)+x(4)^2))-2*beta*x(1)*x(2)*x(5))*(1-beta*x(1)^2/2)
    x(7)
    -Pr*x(1)*x(7)
    x(9)
    Sc*(sigma*(1+delta*x(6))^n*exp(-E/(1+delta*x(6)))*x(8)-x(1)*x(9))];
end

function [t,x]=euler(fun, x0, tend, dt)
t=0:dt:tend;
n=tend/dt;
x=x0;
for i=1:n
    x(:, i+1)=x(:,i)+fun(t,x(:,i))*dt;
end
end