 clc; clear all;
% close all;
%% BARRIER Input parameters%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B=0.0001; %Mainland slope
B=10^100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dt=10;% Toe depth (meters). Typically in the range 10-20m
We=800; %Equilibrium width (meters)
He=2;  %Equilibrium heigth (meters)
Ae=0.02; %Equilibrium shoreface slope
Qow_max=10; %Maximum overwash flux (m^2/year)
Vd_max=We*He; %Maximum deficit volume (m^2/year)
K=2000;  %Shoreface Flux constant (m^2/year)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zdot=0.005; %sea-level rise rate [m/yr]
bm1c=1000; %Critical marsh width (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MARSH Input parameters%%%%%%%%%%%%%%%
rhos=1000;rhom=1000; %marsh and mudlfat density. must be the same to conserve mass
P=12.5/(24*365); %tidal period, yr
ws=0.5*10^-3*(365*24*3600); %settling velocity, m/yr  (m/yr, multiplied my the conversion factor...)
lamda=0.0001; %parameter for sed. conc., adim
dist=10; %parameter for marsh erosion, m
amp=1.4/2; %tidal amplitude, equal to half the tidal range, m
% amp=2/2; %tidal amplitude, equal to half the tidal range, m
rng=2*amp;
Dmin=0;
Dmax=.7167*rng-.0483;
% delT=0;sigB=0.06;%degrees C-1=
tcr=0.1;  %Pa
wind=10; %m/s
Ka=2; %marsh progradation coeff. adimensional
% Ke=0.15; %bank erosion coeff. 
Ke=0.15; %bank erosion coeff. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Co=0.1;  %external sed. conc. [g/l =kg/m3]
%%% SALT MARSH GROWTH%%%%%%%%%%%%%%%%%%%%%%%%%%%5
BMax=2.500;% kg/m2
AA=.25.*(Dmax-Dmin)^2;por=1000/2650;
chiref=0.158;% Only the refractory part
po=1000; %dens of org
%% BARRIER Initial conditions %%%%%%%%%%%%
A=1.0*Ae;W=1.0*We;H=1.0*He; %Barrier initially in equlibrium
Ao=A;Wo=W;Ho=H; 
xt=0;xs=Dt/A;xb=xs+W;
xto=xt;xso=xs;xbo=xb;
Z=1.0*Dt;Zo=1.0*Dt;
%% Initial marsh-lagoon configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bfo=10000;
bm1o=1000;
bm2o=2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BL=bfo+bm1o+bm2o;
xm1o=xbo+bm1o;
xm2o=xbo+bm1o+bfo;
xL=xbo+bm1o+bfo+bm2o;
dfo=2;% Depth lagoon respect to MHW
dmo=Dmax/2;% Depth marsh respect to MHW
bf=bfo; 
bm1=bm1o;bm2=bm2o; 
xm1=xm1o; xm2=xm2o; 
df=dfo;dm1=dmo;dm2=dmo;dm=dm1;
xm3o=xm2o+bm2o-(dm-amp)/B;
xm3=xm3o;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computational parameters%%%%%%%%%%%%%%%
Tmax=1000; %Runing time (years)
dt=0.01;    
t=0:dt:Tmax;
n=length(t);
Heigth=zeros(1,n); Width=zeros(1,n);OverwashFlux=zeros(1,n); 
XB=zeros(1,n);XT= zeros(1,n);XS=zeros(1,n);
Heigth(1)=He; Width(1)=He;XB(1)=xb;XT(1)=xt;XS(1)=xs;
SHrate=zeros(1,n);XTrate=zeros(1,n);XBrate=zeros(1,n);
XM1=zeros(1,n);XM2=zeros(1,n);DF=zeros(1,n);
DM1=zeros(1,n);DM2=zeros(1,n);
BM1=zeros(1,n);BM2=zeros(1,n);
XL=zeros(1,n);BF=zeros(1,n);ZZ=zeros(1,n);
QOWBM=zeros(1,n);QOWBL=zeros(1,n);QOWB=zeros(1,n);
WO=zeros(1,n);
XBDOT=zeros(1,n);
DMDOT=zeros(1,n);
DDD=zeros(1,n);
XM3=zeros(1,n);
CR=zeros(1,n);
tdrown_W=0;
t_Lfill=0;
%% Main code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
Z=Z+zdot*dt;
%% excess width. Note: if bm>DW, then no sediment goes on the marsh progradation. 
Phi=min(1,bm1/bm1c); %relative accomodation width. How much goes on the marsh
EDm=dm1-amp; %top of marsh
EDl=df-amp; %lagoon
DD=Phi*EDm+EDl*(1-Phi);
%% Deficit Volume
Vd_H=max((He-H)*W,0);
Vd_B=max((We-W)*(H+DD),0);
Vd=Vd_H+Vd_B;
%% Overwash formulation
if Vd<Vd_max; Qow_H=Qow_max*Vd_H/Vd_max;Qow_B=Qow_max*Vd_B/Vd_max;else
Qow_H=Qow_max*Vd_H/Vd;Qow_B=Qow_max*Vd_B/Vd;end
Qow=Qow_H+Qow_B;
Qsf=K*(Ae-A);
Qow_B_m=Qow_B*Phi;
Qow_B_l=Qow_B*(1-Phi);
%% BARRIER EQUATIONS
Hdot=Qow_H/W-zdot; 
xbdot=Qow_B_m/(H+EDm);
xsdot=2*Qow/(Dt+2*H)-4*Qsf*(H+Dt)/(2*H+Dt)^2;
xtdot=2*Qsf*(1/(Dt+2*H)+1/Dt)+2*zdot/A;
H=H+Hdot*dt;if H<0;tdrown_H = t(i);jj=i-1; break; end
xb=xb+xbdot*dt;
xs=xs+xsdot*dt;
xt=xt+xtdot*dt;
A=Dt/(xs-xt);
W=xb-xs;if W<0;tdrown_W = t(i); jj=i-1; break; end
% Variable storage %%%%%%%%%%%%%%%%%%%%
Heigth(i)=H; Width(i)=W; OverwashFlux(i)=Qow;
XB(i)=xb;XT(i)=xt;XS(i)=xs;ZZ(i)=Z;
QOWBM(i)=Qow_B_m;QOWBL(i)=Qow_B_l;QOWB(i)=Qow_B;
WO(i)=Wo;DDD(i)=DD;
XBrate(i)=xbdot;

%% BACKBARRIER EQUATIONS
%% salt marsh organic deposition %%%%%%%%%%%%%%%%%%%%%%%
Bpeak=BMax*(Dmax-dm)*(dm-Dmin)/AA;if (Bpeak<=1e-3);Bpeak=0; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bfrac=(Bpeak/BMax);
nuGp=.0138;
AMC=(180.)*Bpeak*(nuGp);
Rref=AMC*chiref;
FFm=(1/por)*(Rref/po);
%%% average depths
Df=(df+(df-min(rng,df)))/2; 
Dm=(dm+(dm-min(rng,dm)))/2; 
%%% evolution of marsh and mudflat
tw= wavetau(bf,wind,Df);
if Dm>1e-4;twm=wavetauBmod(bf,wind,Dm ,Bfrac);else twm=0;end
tau=max((tw-tcr)/tcr,0)*lamda;
taum=max((twm-tcr)/tcr,0)*lamda;
Cr=rhos*tau/(1+tau);
Cm=rhos*taum/(1+taum);
hb=dm+(df-dm)*(1-exp(-dist*0.1/df)); %scarp height according to profile
WP=waveTRNS(amp,wind,bf,hb);  %wave power
Fc=(Cr-Co)*min(rng,df)/P/rhom;%flux mudflat - continental shelf (mud)
E=(Ke*WP/(hb-dm)-Ka*Cr*ws/rhom);
Fm=(Cr-Cm)*min(rng,dm)/P/rhom; %flux from tidal flat to marsh
%% Moving Boundaries
%%% Rates
dmdot=-Fm-FFm+zdot;
xm1dot=-E+Qow_B_l/(df-dm); 
xm2dot=E;
xm3dot=(zdot-dmdot)/B;
dfdot=-E*(df-dm1)/bf-E*(df-dm2)/bf+Fm*(bm1+bm2)/bf+Fc+zdot;
%%% Absolute value
xm1=xm1+xm1dot*dt;
xm2=xm2+xm2dot*dt;
% xm3=xm2o+bm2o+(Z-Zo-(dm-amp))/B;
xm3=xm3+xm3dot*dt;
dm=dm+dmdot*dt;
df=df+dfdot*dt;
xL=xm2o+bm2o+(Z-Zo)/B; %% xL is equivalent to xm3.
%% Marsh 1 or marsh 2 eroded away %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if xm1<=xb;xm1=xb; end
if xm3<=xm2;xm2=xm3; end
%% Filling of the lagoon or Vertical marsh collapse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dm1>Dmax; t_vmc=t(i);break; end %% I need to redistribute the sediments 
if df<=dm1; t_Lfill=t(i);xm1=xb; xm2=xm3; end
dm1=dm;
dm2=dm;
bm1=xm1-xb;
bm2=xm3-xm2;
bf=xm2-xm1;
% if bf>BL; break; end
%% Variable storage %%%%%%%%%%%%%%%%%%%%
XM1(i)=xm1;XM2(i)=xm2;XM3(i)=xm3; 
DF(i)=df;DM1(i)=dm1;DM2(i)=dm2;BF(i)=bf;
BM1(i)=bm1;BM2(i)=bm2;XL(i)=xL;
CR(i)=Cr;
DMDOT(i)=dmdot;
end

%% SUBPLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
hold on
subplot(2,2,1);
hold on
plot(t,BM1,'k', 'linewidth',2)
% plot(t,BM2,'k--', 'linewidth',2)
box on
xlabel('time [years]')
title('Backbarrier marsh width (m)')

subplot(2,2,3);
hold on
plot(t,Width,'k', 'linewidth',2)
xlabel('time [years]')
if tdrown_W>0; xlim([0 tdrown_W]); end
title('Barrier Width (m)')
box on

subplot(2,2,2);
hold on
plot(t,BF,'k', 'linewidth',2)
xlabel('time [years]')
if tdrown_W>0; xlim([0 tdrown_W]); end
title('Lagoon width (m)')
box on

subplot(2,2,4);
hold on
plot(t,DF,'k', 'linewidth',2)
xlabel('time [years]')
if tdrown_W>0; xlim([0 tdrown_W]); end
title('Lagoon depth (m)')
box on




