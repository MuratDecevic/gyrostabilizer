
function eqns = noGyrodynamics(t,x,p)
% diff equations describing the dynamics of the vessel

Isys=p(1); B44=p(2); C44=p(3);


th=x(1);
thp=x(2);
%F40=400*sin(2*pi*f*t);
F40=waveForcing(t);
% equations
%Isys*theta''+B44*theta'+C44*theta=F40

thpp=(F40/Isys)-(B44/Isys)*thp-(C44/Isys)*th;


eqns=[thp;thpp];