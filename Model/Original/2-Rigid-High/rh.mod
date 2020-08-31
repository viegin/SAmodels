% Monetary Policy In An Economy With High Structural Unemployment
% DSGE Model

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Definition of variables
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

var pi, ahat, uhat, chat, i, dm, m, nhat, xhat, l, w, h; 
varexo ea, ed, em, el;

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
parameters alp, bet, phi, chi, lam, gam, theta, e, B, M, x, g, u, delt, rhoa, rhod, rhol, rhom, psi, rho, net0, net1, h0, hL, hF, k0, kL, kF, phip, phic, phiu;
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% parameter values:

alp = 1; 						% a constant -- elasticity of hiring cost
bet = 0.99; 					% discount rate 
theta=0.75;
gam = 0.5; 						% index of wage rigidity
e = 6; 							% elasticity of subsititution
B = 0.12; 						% level of hiring cost
M = 1.2; 						% optimal gross markup
x = 0.15; 						% labour market tightness index
g = B*x^alp; 					% hiring cost
u = 0.3; 						% unemployment
phi=1;
chi=(1-(1-bet*(1-delt))*(1+alp)*B*x^alp-bet*(1-delt)*alp*B*x^(1+alp))/((1-delt*B*x^alp)*(x/(delt+(1-delt)*x)));
delt = u*x/((1-u)*(1-x)); 		% exogenous separation rate
psi = 1-(1-bet*(1-delt))*g*M; 	% Philips curve parameter
rhoa = 0.9; 					% persistence of the productivity process
rhod = 0.9;
rhol = 0.9;
rhom = 0.9;
rho = -log(bet);
lam =((1-bet*theta)*(1-theta))/theta; 					% price duration
net0 = (1-g*(1+alp))/(1-delt*g); 
net1 = g*(1-delt)*(1+alp*(1-x))/(1-delt*g);
h0 = (alp*g*M/delt)*(1+bet*(1-delt)*(1-delt)*(1-x))+bet*(1-delt)*g*M*(net1-net0);
hL = -(alp*g*M/delt)*(1-delt)*(1-x)-bet*(1-delt)*g*M*net1; 
hF = -bet*(1-delt)*g*M*((alp/delt)-net0);
k0 = lam*h0/(1-u); 
kL = -lam*hL/(1-u);
kF = -lam*hF/(1-u); 
phip = 1.5; 
phic = 0.5; 
phiu = 0; 

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Model
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

model;

%uhat = 1-nhat(-1)+delt*nhat(-1);
h=nhat-nhat(-1)-delt*nhat(-1);
%xhat=h-uhat;
%uhat=1-nhat;
log((1-(1-bet*(1-delt))*(1+alp)*B*x^alp-bet*(1-delt)*alp*B*x^(1+alp))/((1-delt*B*x^alp)*(x/(delt+(1-delt)*x))))+phi*nhat=w-pi-chat;


pi = bet*pi(+1)-k0*uhat+kL*uhat(-1)+kF*uhat(+1)-lam*psi*gam*ahat; 										% Phillips curve
delt*xhat = nhat - (1-delt)*(1-x)*nhat(-1); 															% relation between labour market tightness and employment
chat = ahat +((1-g)/(1-delt*g))*nhat + (g*(1-delt)/(1-delt*g))*nhat(-1) - (alp*g/(1-delt*g))*delt*xhat+l; % consumption
chat = chat(+1) - (i - pi(+1) - rho)+dm; 																% consumer's first order condition
uhat=-(1-u)*nhat; 																						% relation between unemployment and employment
i = phip*pi+phic*chat+phiu*uhat+m; 																		% taylor rule

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% exogenous process:
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ahat = rhoa*ahat(-1)+ea;
dm = rhod*dm(-1)+ed;
m = rhom*m(-1)+em;
l = rhol*l(-1)+el;


end;


initval;
pi = 0;
uhat = 0;
ahat = 0;
chat = 0;
i = 0;
dm = 0;

end;


steady;
check;

shocks;
var ea; stderr 1;
var ed; stderr 1;
var em; stderr 1;
var el; stderr 1;
end;


stoch_simul(TeX);
%collect_latex_files;

%if system(['pdflatex -halt-on-error ' M_.fname '_TeX_binder.TeX'])
%     error('TeX-File did not compile.')
% end
