clearvars
clear all

%defining constants
N = 100;
R = 8.314;
P1 = 10;
T = 298;
L = 1e-4;
dx = L/N;
dt = 1e-7;
D_ab = 1e-6;
x = linspace(0,L,N);
Ci = zeros(N,1);
Cf = zeros(N,1);
Ci(1) = 0;
C = zeros(N,1);
Co = P1/(R*T);
C(1) = Co;
C(N) = 0;
t_max = [1e-5,0.0001,0.0005,0.001,0.01];
for k = 1:5
    for i = 0:dt:t_max(k)
        Cf(1) = P1/(R*T);
        for j = 2:N-1
            Cf(j) = Ci(j) +     ...
            ((dt*D_ab)/(dx)^2)*(Ci(j+1) ...
            - 2*Ci(j) + Ci(j-1));
            temp_var = x(j) / (2*sqrt(D_ab*i));
            C(j) = Co*(1 - erf(temp_var));
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
legend('1e-5s','1e-5s(analytical)','1e-4s','1e-4s(analytical)','5e-4s','5e-4s(analytical)','0.001s','0.001s(analytical)','0.01s','0.01s(analytical)')