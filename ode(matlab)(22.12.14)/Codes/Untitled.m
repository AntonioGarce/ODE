% Plot Figure 2
clc; clear; close all

global beta lambda Pr sigma n Sc E delta x0 
beta =0.5;
lambda=0.2; 
Pr = 7;
sigma = 1;
n = 0.5;
Sc =1;
E =5;
delta = 1;
etaend = 4;

% u(1) = -1.1;
% u(2) = 0.5;
% u(3) = -0.51903;
% u(4) = -0.95064;
u(1) = -1.13;
u(2) = -0.23;
u(3) =-1.84545;
u(4) = -0.55;

% u=fsolve(@Fshoot, u)
% Fshoot(u)

x0 = [0; 1; u(1); 0 ; u(2); 1; u(3); 1; u(4)];
solinit=bvpinit(linspace(0, 12, 101), @mat4init);
sol = bvp4c(@odefun20, @mat4bc, solinit);
x=linspace(0, 12, 101);
y=deval(sol,x);
plot(x,y(4,:))

function xinit = mat4init(eta)
global x0;
xinit = [x0(1)+x0(2)*eta%+x0(3)*eta^2/2
    x0(2)+x0(3)*eta
    x0(3)
    x0(4)+x0(5)*eta
    x0(5)
    x0(6)+x0(7)*eta
    x0(7)
    x0(8)+x0(9)*eta
    x0(9)
    ];
end
