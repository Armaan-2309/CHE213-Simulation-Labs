clc
clear

R = 0.015;
tmax = 100*60;
Dp = 0.000000195;
taumax = Dp*tmax/(R^2);
tspan = [0,taumax];
y0 = [0.99,0.99];
[t,y] = ode45(@(t,y) odefun(t,y), tspan, y0);
Co = 0.02115;
Ct = y(:,2)*Co;
Rf = y(:,1)*R;
time = t*(R^2)/Dp;

figure
plot(time,Ct);
xlabel('Time(in s)');
ylabel('Concentration');
title('Concentration as a function of time');
figure
plot(time,Rf);
xlabel('Time');
ylabel('Radius');
title('Radius as a function of time');


%%On Increasing Stirring Speed


R = 0.015;
tmax = 100*60;
Dp = 0.000000195;
taumax = Dp*tmax/(R^2);
tspan = [0,taumax];
y0 = [0.99,0.99];
[t,y] = ode45(@(t,y) odefun(t,y), tspan, y0);
Co = 0.02115;
Ct = y(:,2)*Co;
Rf = y(:,1)*R;
time = t*(R^2)/Dp;
y0 = [0.99,0.99];
[t,y] = ode45(@(t,y) odefunc_a(t,y), tspan, y0);
Ct1 = y(:,2)*Co;
Rf1 = y(:,1)*R;
time_a = t*(R^2)/Dp;

figure
plot(time,Ct);
hold on
plot(time_a,Ct1);
xlabel('Time(in s)');
ylabel('Concentration');
title('Reduced stirring speed');
legend("k=6.6276e-005","k=6.6276e-004")


%%On reducing initial Adsorbent amount


R=0.015;
tmax=100*60;
Dp=0.000000195;
taumax = Dp*tmax/(R^2);
tspan = [0,taumax];
y0 = [0.99,0.99];
[t,y] = ode45(@(t,y) odefun(t,y), tspan, y0);
Co=0.02115;
Ct = y(:,2)*Co;
Rf = y(:,1)*R;
time = t*(R^2)/Dp;
y0 = [0.99,0.99];
[t,y] = ode45(@(t,y) odefunc_b(t,y), tspan, y0);
Ct2 = y(:,2)*Co;
Rf2 = y(:,1)*R;
time_b = t*(R^2)/Dp;

figure
plot(time,Ct);
hold on
plot(time_b,Ct2);
xlabel('Time(in s)');
ylabel('Concentration');
title('Reduced initial adsorbent amount');
legend("W=8","W=7")


%%On Increasing Effective Pore Diffusivity


R=0.015;
tmax=100*60;
Dp=0.000000195;
taumax = Dp*tmax/(R^2);
tspan = [0,taumax];
y0 = [0.99,0.99];
[t,y] = ode45(@(t,y) odefun(t,y), tspan, y0);
Co=0.02115;
Ct = y(:,2)*Co;
Rf = y(:,1)*R;
time = t*(R^2)/Dp;
y0 = [0.99,0.99];
[t,y] = ode45(@(t,y) odefunc_c(t,y), tspan, y0);
Ct3 = y(:,2)*Co;
Rf3 = y(:,1)*R;
time_c = t*(R^2)/Dp;

figure
plot(time,Ct);
hold on
plot(time_c,Ct3);
xlabel('Time(in s)');
ylabel('Concentration');
title('Higher effective pore diffusivity');
legend("Dp=0.000000195*2","Dp=0.000000195")


%Question 2(a)

R = 0.015;
tmax = 200*60;
Dp = 0.000000195;
taumax = Dp*tmax/(R^2);
tspan = [0,taumax];
y0 = [0.99,0.99];
[t,y] = ode45(@(t,y) odefun(t,y), tspan, y0);
Co = 0.02115;
Ct = y(:,2)*Co;
Rf = y(:,1)*R;
time = t*(R^2)/Dp;
y0 = [0.99,0.99];
[t,y] = ode45(@(t,y) odefunc_a(t,y), tspan, y0);
Ct1 = y(:,2)*Co;
Rf1 = y(:,1)*R;
time_a = t*(R^2)/Dp;

figure
plot(time,Ct);
hold on
plot(time_a,Ct1);
xlabel('Time(in s)');
ylabel('Concentration');
title('Increased stirring speed at eqb');
legend("k=6.6276e-005","k=6.6276e-004")


%Question 2(b)


R=0.015;
tmax= 200*60;
Dp=0.000000195;
taumax = Dp*tmax/(R^2);
tspan = [0,taumax];
y0 = [0.99,0.99];
[t,y] = ode45(@(t,y) odefun(t,y), tspan, y0);
Co=0.02115;
Ct = y(:,2)*Co;
Rf = y(:,1)*R;
time = t*(R^2)/Dp;
y0 = [0.99,0.99];
[t,y] = ode45(@(t,y) odefunc_c(t,y), tspan, y0);
Ct3 = y(:,2)*Co;
Rf3 = y(:,1)*R;
time_c = t*(R^2)/Dp;

figure
plot(time,Ct);
hold on
plot(time_c,Ct3);
xlabel('Time(in s)');
ylabel('Concentration');
title('Higher effective pore diffusivity at eqb');
legend("Dp=0.000000195*2","Dp=0.000000195")

%%%
%%Function:
function dydt = odefun(t,y)
 kf=6.6276e-005; %liquid phase mass transfer coefficient(cm/s)
 R=0.015; %Adsorbent particle radius(cm)
 Dp=0.000000195; %Effective diffusion coefficient in the adsorbent(cm2/s)
 Co=0.02115; %initial liquid phase concentration(g/L)
 V=0.4; %Volume of batch reactor(L)
 W=8.0; %weight of the adsorbent(g)
 ro=1350; %adsorbent density
 Ch=W/(V*Co);
 Bi=(kf*R)/Dp;
 ys=0.580*2.90 ; %Langmuir isotherm constant(L/g)
 ko=2900 ; %Langmuir isotherm constant(L/g)
 ko1=Co*ko; %From Non-dimensional Langmuir eqn
 yes=Co*ys; %From Non-dimensional Langmuir eqn
 tmax=6000; %Experiment

 Cet_d = @(Ct_d,r) Bi*(1-r)*Ct_d/(r + Bi*(1-r));

 Yet = @(Ct_d,r) yes*Cet_d(Ct_d,r)/(1 + ko1*Cet_d(Ct_d,r));

 f1 = @(Ct_d,r) -1*Bi*(Co/(ro*Yet(Ct_d,r)))*(Ct_d - Cet_d(Ct_d,r))/(r^2);

 M = @(Ct_d,r) 1 + Ch*(1 - r^3)*yes*Bi*(1-r)/(((1 + ko1*Cet_d(Ct_d,r))^2)*(r + (1-r)*Bi));

 N = @(Ct_d,r) 3*Ch*Yet(Ct_d,r)*(r^2) + Ch*yes*Bi*(1 - r^3)*Ct_d/(((1 + ko1*Cet_d(Ct_d,r))^2)*((r + (1-r)*Bi)^2));
 %y(1) = r; y(2) = Ct_d;
 dydt = [f1(y(2),y(1)); N(y(2),y(1))*f1(y(2),y(1))/M(y(2),y(1))];
end
function dydt = odefunc_a(t,y)
 kf=6.6276e-004; %liquid phase mass transfer coefficient(cm/s)
 R=0.015; %Adsorbent particle radius(cm)
 Dp=0.000000195; %Effective diffusion coefficient in the adsorbent(cm2/s)
 Co=0.02115; %initial liquid phase concentration(g/L)
 V=0.4; %Volume of batch reactor(L)
 W=8.0; %weight of the adsorbent(g)
 ro=1350; %adsorbent density
 Ch=W/(V*Co);
 Bi=(kf*R)/Dp;
 ys=0.580*2.90 ; %Langmuir isotherm constant(L/g)
 ko=2900 ; %Langmuir isotherm constant(L/g)
 ko1=Co*ko; %From Non-dimensional Langmuir eqn
 yes=Co*ys; %From Non-dimensional Langmuir eqn
 tmax=6000; %Experiment

 Cet_d = @(Ct_d,r) Bi*(1-r)*Ct_d/(r + Bi*(1-r));

 Yet = @(Ct_d,r) yes*Cet_d(Ct_d,r)/(1 + ko1*Cet_d(Ct_d,r));

 f1 = @(Ct_d,r) -1*Bi*(Co/(ro*Yet(Ct_d,r)))*(Ct_d - Cet_d(Ct_d,r))/(r^2);

 M = @(Ct_d,r) 1 + Ch*(1 - r^3)*yes*Bi*(1-r)/(((1 + ko1*Cet_d(Ct_d,r))^2)*(r + (1-r)*Bi));

 N = @(Ct_d,r) 3*Ch*Yet(Ct_d,r)*(r^2) + Ch*yes*Bi*(1 - r^3)*Ct_d/(((1 + ko1*Cet_d(Ct_d,r))^2)*((r + (1-r)*Bi)^2));
 %y(1) = r; y(2) = Ct_d;
 dydt = [f1(y(2),y(1)); N(y(2),y(1))*f1(y(2),y(1))/M(y(2),y(1))];
end


function dydt = odefunc_b(t,y)
 kf=6.6276e-005; %liquid phase mass transfer coefficient(cm/s)
 R=0.015; %Adsorbent particle radius(cm)
 Dp=0.000000195; %Effective diffusion coefficient in the adsorbent(cm2/s)
 Co=0.02115; %initial liquid phase concentration(g/L)
 V=0.4; %Volume of batch reactor(L)
 W=7.0; %weight of the adsorbent(g)
 ro=1350; %adsorbent density
 Ch=W/(V*Co);
 Bi=(kf*R)/Dp;
 ys=0.580*2.90 ; %Langmuir isotherm constant(L/g)
 ko=2900 ; %Langmuir isotherm constant(L/g)
 ko1=Co*ko; %From Non-dimensional Langmuir eqn
 yes=Co*ys; %From Non-dimensional Langmuir eqn
 tmax=6000; %Experiment

 Cet_d = @(Ct_d,r) Bi*(1-r)*Ct_d/(r + Bi*(1-r));

 Yet = @(Ct_d,r) yes*Cet_d(Ct_d,r)/(1 + ko1*Cet_d(Ct_d,r));

 f1 = @(Ct_d,r) -1*Bi*(Co/(ro*Yet(Ct_d,r)))*(Ct_d - Cet_d(Ct_d,r))/(r^2);

 M = @(Ct_d,r) 1 + Ch*(1 - r^3)*yes*Bi*(1-r)/(((1 + ko1*Cet_d(Ct_d,r))^2)*(r + (1-r)*Bi));

 N = @(Ct_d,r) 3*Ch*Yet(Ct_d,r)*(r^2) + Ch*yes*Bi*(1 - r^3)*Ct_d/(((1 + ko1*Cet_d(Ct_d,r))^2)*((r + (1-r)*Bi)^2));
 %y(1) = r; y(2) = Ct_d;
 dydt = [f1(y(2),y(1)); N(y(2),y(1))*f1(y(2),y(1))/M(y(2),y(1))];
end


function dydt = odefunc_c(t,y)
 kf=6.6276e-005; %liquid phase mass transfer coefficient(cm/s)
 R=0.015; %Adsorbent particle radius(cm)
 Dp=0.000000195*2; %Effective diffusion coefficient in the adsorbent(cm2/s)
 Co=0.02115; %initial liquid phase concentration(g/L)
 V=0.4; %Volume of batch reactor(L)
 W=8.0; %weight of the adsorbent(g)
 ro=1350; %adsorbent density
 Ch=W/(V*Co);
 Bi=(kf*R)/Dp;
 ys=0.580*2.90 ; %Langmuir isotherm constant(L/g)
 ko=2900 ; %Langmuir isotherm constant(L/g)
 ko1=Co*ko; %From Non-dimensional Langmuir eqn
 yes=Co*ys; %From Non-dimensional Langmuir eqn
 tmax=6000; %Experiment

 Cet_d = @(Ct_d,r) Bi*(1-r)*Ct_d/(r + Bi*(1-r));

 Yet = @(Ct_d,r) yes*Cet_d(Ct_d,r)/(1 + ko1*Cet_d(Ct_d,r));

 f1 = @(Ct_d,r) -1*Bi*(Co/(ro*Yet(Ct_d,r)))*(Ct_d - Cet_d(Ct_d,r))/(r^2);

 M = @(Ct_d,r) 1 + Ch*(1 - r^3)*yes*Bi*(1-r)/(((1 + ko1*Cet_d(Ct_d,r))^2)*(r + (1-r)*Bi));

 N = @(Ct_d,r) 3*Ch*Yet(Ct_d,r)*(r^2) + Ch*yes*Bi*(1 - r^3)*Ct_d/(((1 + ko1*Cet_d(Ct_d,r))^2)*((r + (1-r)*Bi)^2));
 %y(1) = r; y(2) = Ct_d;
 dydt = [f1(y(2),y(1)); N(y(2),y(1))*f1(y(2),y(1))/M(y(2),y(1))];
end