%%

% inputs
V = 0.05; % voltage in Volts
s = 0.3; % how many electronic charges per base
d = 0.13; % diffusion constant in units of D_0
T = 273.15 + 25; % K
n = 14; % number of bases past constriction
tolerance = 1e-3; % for numerical integration, accurate to this

% constants and derived quantities
e_0 = 1.6e-19; % coulombs
nu = 0.59; % Flory exponent
D_0 = 3.31e-10; % diffusion const of 1.5nm sphere, m^2/s
D = d * D_0;
kB = 1.38e-23; % SI units
Lp = 1e-9; % Tinland 1997
Lb = 0.5e-9; % length per base ssDNA, meters
L = n * Lb; % length of domain, meters
sigma = s * e_0 / Lb; % charge per unit length

f = @(x,v) exp(sigma*v*x/(kB*T));
xx = linspace(0,L,1000);

figure
plot(xx*1e9*2,f(xx,-0.08))
hold on
plot(xx*1e9*2,f(xx,-0.16))
plot([2,4],f(4*Lb,-0.08)*ones(1,2),'o','Color',[0.5,0.5,0.5])
line([2,4],f(4*Lb,-0.08)*ones(1,2),'LineStyle','--','Color',[0.5,0.5,0.5])
plot([3,6],f(6*Lb,-0.08)*ones(1,2),'o','Color',[0.5,0.5,0.5])
line([3,6],f(6*Lb,-0.08)*ones(1,2),'LineStyle','--','Color',[0.5,0.5,0.5])
plot([4,8],f(8*Lb,-0.08)*ones(1,2),'o','Color',[0.5,0.5,0.5])
line([4,8],f(8*Lb,-0.08)*ones(1,2),'LineStyle','--','Color',[0.5,0.5,0.5])
text(1,1e-3,'160mV')
text(5,1e-3,'80mV')
text(2,1e-2,'2 bases')
text(3,1e-4,'3 bases')
text(6,1e-5,'4 bases')
set(gca,'yscale','log')
set(gca,'fontsize',12,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])
ylabel('Fraction of time spent past given base')
xlabel('Base number')
ylim([1e-7 1])
xlim([0 10])
title([num2str(s) 'e^- per base'])
set(gcf,'position',[-585   463   478   490])
