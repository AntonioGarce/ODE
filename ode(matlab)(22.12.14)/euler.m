function [t,x]=euler(fun, x0, tend, dt)
t=0:dt:tend;
n=tend/dt;
x=x0;
for i=1:n
    x(:, i+1)=x(:,i)+fun(t,x(:,i))*dt;
end

