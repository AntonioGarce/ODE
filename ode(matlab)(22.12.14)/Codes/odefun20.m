function dxdeta = odefun20(eta,x)
global beta lambda Pr sigma n Sc E delta
dxdeta = [x(2);
    x(3);
    (x(2)^2-x(1)*x(3)-2*lambda*(x(4)-beta*x(1)*x(5))-2*beta*x(1)*x(2)*x(3))/(1-beta*x(1)^2/2);
    x(5);
    (x(2)*x(4)-x(1)*x(5)+2*lambda*(x(2)+beta*(x(2)^2-x(1)*x(3)+x(4)^2))-2*beta*x(1)*x(2)*x(5))*(1-beta*x(1)^2/2);...
    x(7);
    -Pr*x(1)*x(7);
    x(9);
    Sc*(sigma*(1+delta*x(6))^n*exp(-E/(1+delta*x(6)))*x(8)-x(1)*x(9))];
end

