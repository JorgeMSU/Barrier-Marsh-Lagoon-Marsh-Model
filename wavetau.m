function tw=wavetau(fetch,wind,Df)
[Hs,Tp]=YeV(fetch,wind,Df);
kk=wavek(1./Tp,Df);
Um=(pi*Hs./Tp./sinh(kk.*Df));
aw=Tp*Um/pi;ko=0.001;
fw=0.4*(aw/ko)^-0.75;
tw=1/2*1020*fw*Um^2;