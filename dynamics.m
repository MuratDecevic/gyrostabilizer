
function eqns = dynamics(t,x,p)
% diff equations describing the dynamics of the vessel and gyroscope

Isys=p(1); B44=p(2); C44=p(3);
Is=p(4); ws=p(5); It=p(6); Bg=p(7);
Cg=p(8);

th=x(1); % roll angle, rad
thp=x(2); % roll rate, rad/s
b=x(3); % precession angle, rad
bp=x(4); % precession rate, rad/s

F40=waveForcing(t); % wave forcinf term
% equations
%It*beta''+Bg*beta'+Cg*sin(beta)=Is*ws*theta'*cos(beta)+(It/Is)*beta'
%Isys*theta''+B44*theta'+C44*theta=F40-Is*ws*beta'*cos(beta)
% The beta'/Is term is the precession control input.
thpp=(F40/Isys)-(Is/Isys)*ws*bp*cos(b)-(B44/Isys)*thp-(C44/Isys)*th;
bpp=(Is/It)*ws*thp*cos(b)-(Bg/It)*bp-(Cg/It)*sin(b)+bp/Is;

eqns=[thp;thpp;bp;bpp];
