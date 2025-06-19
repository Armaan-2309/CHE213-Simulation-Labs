
clc
clear all


D_ab=1e-6;
P1=10;
R=8.314;
T=298;
N=100;
L=1e-4;
dt=1e-8;
dx=L/N;

x=linspace(0,L,N);
C_old=zeros(N,1);
C_new=zeros(N,1);
C_old(1)=P1/R*T;
C_old(N)=0;
C_new(1)=P1/R*T;
C_old(N)=0;
C_old1=zeros(N,1);
C_new1=zeros(N,1);
C_old1(1)=P1/R*T;
C_old1(N)=0;
C_new1(1)=P1/R*T;
C_old1(N)=0;
u=0.1;

C_old(1)=P1/(R*T);

times=[1e-5,1e-4,1e-3];
for k=1:3
    for i=0:dt:times(k)
        for j=2:N-1
            C_new(j)=C_old(j)-(((u*dt)/(dx))*(C_old(j)-C_old(j-1)))+(dt/(dx).^2)*(D_ab)*(C_old(j+1)-2*C_old(j)+C_old(j-1));%iterations for first method
            C_new1(j)=C_old1(j)+(dt/(dx).^2)*(D_ab)*(C_old1(j+1)-2*C_old1(j)+C_old1(j-1));%iterations for first method
        end

        C_old=C_new;
        C_old1=C_new1;
    end
figure(1);
plot(x,C_new);
hold on;
plot(x,C_new1,"o");
hold on;
title("Concentration profiles(diffusion v/s convection");
xlabel("Position(m)");
ylabel("Concentration(mol/m^3)");
legend("t = 1e-5s","t = 1e-5s(wc)","t = 1e-4s","t = 1e-4s(wc)","t = 1e-3s","t = 1e-3s(wc)");
end
%Both methods agree initially but they diverge at later instances of times.