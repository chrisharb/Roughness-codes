function [Z , Lcc, Lct, Lacc, Lact] = crustaleq(H, A, G)

%%
% This is a function to compute seismicity distributions on rough
% continental scale faults, based on the geometric energy balance of
% Harbord, Nielsen & DePaola (2016). 
% Input parameters:
%           H  = The hurst exponent of a self-affine fault surface
%           A  = The amplitude scaling of a self-affine fault surface
%           G  = The surface fracture energy of a typical earthquake, [J
%           m^-2]
%%
rho = 2600; %Average density of the continental cust
alpha = 5800; %P-wave velocity
beta = 3200; %S-wave velocity
g = 9.81; %Acceleration due to gravity
Z = linspace(1000,30000,1000); %Depth model
Sz = rho*g.*Z; %Calculate the vertical component of stress
theta = 60; %fault angle relative to Sigma 1
tau0 = 0.15;
taur = 0.1;
C=0;
Mu = beta^2*rho; %Calculate Lames first parameter
lam = alpha^2*rho-2*Mu; %Calculates Lames second parameter

b = (1+sin(tanh(tau0)))./(1-sin(tanh(tau0))); %Conversion for stress component
a = 2*C*b^(1/2); %Conversion for stress component

S1c = a+(b*Sz); %Calculate S1 in compressional environment
S3c = Sz; %Calculate S3 in compressional environment
S1t = Sz; %Calculate S1 in tensional environment
S3t = (Sz-a)./b; %Calculate S3 in tensional environment

Snc = ((S1c+S3c)./2)+(((S1c+S3c)./2).*cos(2*theta)); %Reverse fault normal stress
Snt = ((S1t+S3t)./2)+(((S1t+S3t)./2).*cos(2*theta)); %Normal fault normal stress

Lcc = (8/pi())*((Mu*(lam+Mu))/(lam+2*Mu)).*(G./((tau0.*Snc-taur.*Snc).^2));
Lct = (8/pi())*((Mu*(lam+Mu))/(lam+2*Mu)).*(G./((tau0.*Snc-taur.*Snt).^2));

Lacc = 2*pi()*(2*pi()*(((2*H)/(10^A))^(1/2))^(1/(H-1)).*((Mu./Snc).^(1/(1-H))));
Lact = 2*pi()*(2*pi()*(((2*H)/(10^A))^(1/2))^(1/(H-1)).*((Mu./Snt).^(1/(1-H))));

figure(1)
subplot(1,2,1)
plot(Lcc,Z./1000,Lacc,Z./1000)
%title(['Compressional earthquake scaling']);
xlabel(['Length [m]']);
ylabel(['Depth [km]']);
set(gca,'xscal','log')
set(gca,'Ydir','reverse')
%set(gca,'yscal','log')
legend('Nucleation length','Void scaling')
grid on
subplot(1,2,2)
plot(Lct,Z./1000,Lact,Z./1000)
%title(['Tensional earthquake scaling']);
xlabel(['Length scaling [m]']);
ylabel(['Depth [km]']);
set(gca,'xscal','log')
set(gca,'Ydir','reverse')
set(gca,'yscal','log')

legend('Nucleation length','Void scaling')
grid on
