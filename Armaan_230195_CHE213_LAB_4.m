clc
clear
%% Obtained slope of linear fit curve from excel

m = 0.57085349;

%% Obtained slope of linear fit curve using polyfit

x = [0	0.0186	0.0476	0.0673	0.0881	0.1102	0.1424	0.1894	0.2069];
y = [0	0.0105	0.0272	0.0375	0.0492	0.0624	0.0809	0.1078	0.1182];
X = x./(1 - x);
Y = y./(1 - y);
m = polyfit(X,Y,1);
Yfit = polyval(m,X);
figure
plot(X,Y,'g',LineWidth=1.5)
hold on
plot(X,Yfit,'o')
title('X-Y Equilibrium linear-fit curve')
xlabel('X  (Liq mole ratio)')
ylabel('Y  (Gas mole ratio)')

%% Part (a)
m = m(1);
y = @(x) m*x;
y1 = 0.15;
G_new = 2000/(0.15*46 + 0.85*44);
Gs = G_new*(1 - y1);
y2 = 0.2*0.15*G_new/(Gs + 0.S2*0.15*G_new);
x2 = 0;
Y1 = y1./(1 - y1);
X2 = x2./(1 - x2);
Y2 = y2./(1 - y2);
X1_min = Y1/m;
Ls_min = Gs*(Y1 -Y2)/(X1_min - X2);
Ls = 1.5*Ls_min;
X1 = (Gs/Ls)*(Y1 - Y2) + X2;
S = 0;
Y_eqbm = @(X) m*X;
m_OpLine = Ls/Gs;
Y_OpLine = m_OpLine*X + Y2;
plot(X,Y_OpLine,'b')
Y_temp = Y1;
X_temp = X1;
tol = 0.0001;
while (Y_temp - Y2 >= tol)
    Y_temp = X_temp*m;
    X_temp = (Y_temp - Y2)/m_OpLine;
    S = S + 1;
end
legend('Equilibrium Points','Equilibrium Curve','Operating Line')

disp('Number of stages required is:')
S
