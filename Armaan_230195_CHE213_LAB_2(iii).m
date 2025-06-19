clc
clear all

%defining constants
N = 100;
R = 8.314;
P1 = 10;
T = 298;
L = 1e-4;
dx = L/N;
dt = 1e-8;
D_ab = 1e-6;
x = linspace(0,L,N);
Ci = zeros(N,1);
Cf = zeros(N,1);
Co = P1/(R*T);
Ci(1) = Co;
Ci(N) = 0;
Cf(1) = Co;
Cf(N) = 0;
C = zeros(N,1);
Ci(1) = Co;
% C(N) = 0;
t_max = 1e-2;
u = [0,0.01,0.1,1];
for k = 1:4
    for i = 0:dt:t_max
        Cf(1) = P1/(R*T);
        for j = 2:N-1
            Cf(j) = Ci(j) - (u(k)*(dt/dx))*(Ci(j) - Ci(j-1)) +     ...
            ((dt*D_ab)/(dx)^2)*(Ci(j+1) ...
            - 2*Ci(j) + Ci(j-1));
            temp_var = exp(u(k)*x(j)/D_ab);
            const = exp(u(k)*L/D_ab);
            C(j) = Co*const/(const - 1) + (Co*temp_var)/(1 - const);
        end
        Ci = Cf;
        
        
    end
    plot(x,Cf)
    hold on
    plot(x,C,'p') 
    hold on
end
xlabel("Position  (in m)")
ylabel("Concentration  (in mol/s)")
title("Concentration Profile")
legend('u = 0','u = 0(analyt)','u = 0.01','u = 0.01(analyt)','u = 0.1','u = 0.1(analyt)','u = 1','u = 1(analyt)')