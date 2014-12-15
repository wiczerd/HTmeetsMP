function [excess,theeconomy] = sol_wcPa_ss(wcPa)
% computes a steady state in the economy
% if only one output is requested, it gives excess (as in excess demand)
% if more are requested it outputs:
% theeconomy{1} : N_a
% theeconomy{2} : u
% theeconomy{3} : Q
% theeconomy{4} : J
% theeconomy{5} : Ve
% theeconomy{6} : Vu


global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau

wc = wcPa(1);
Pa = wcPa(2);

J	= (1-lambda)*(Ym-wc)/(1-beta*(1-lambda));
Q	= (kappa/(Ym-wc + beta*J)/Amf)^(-1/eta);
pQ	= Amf*Q^(1-eta);

% solve for Ve Vu as a linear system
flowVe	= (1-lambda)*alpha^alpha*(Pa/(1-alpha)).^(alpha-1).*(wc*(1-tau) -cbar-Pa*abar);
flowVu	= pQ*alpha^alpha*(Pa/(1-alpha)).^(alpha-1).*(wc*(1-tau)-cbar-Pa*abar) ...
		+ (1-pQ)*alpha^alpha*(Pa/(1-alpha)).^(alpha-1).*(be*wc-cbar-Pa*abar);
VeVucoe	= [(1-lambda)*beta, lambda;...
	pQ*beta, (1-pQ)*beta];
VeVucon	= [flowVe;flowVu];

VeVu	= (eye(2)- VeVucoe)\VeVucon;

% now that I have Ve, Vu I can solve for the wage from Nash Bargaing
wc_implied = ( alpha^alpha*(Pa/(1-alpha)).^(alpha-1)*...
	( (1-theta)*be*wc + (theta-1)*beta*(VeVu(1) - VeVu(2)) ) ...
	+ theta*(Ym + beta*J) )...
	/(theta+alpha^alpha*(Pa/(1-alpha)).^(alpha-1)*(1-theta)); 

excess(1) = wc_implied - wc;

% solve for the size of the agricultural economy
% this sets the first order condition of a migrator to be indifferent
Na_supplied = ( ( (1-beta)*VeVu(2)/(alpha^alpha*(Pa/(1-alpha))^(alpha-1)) + cbar +Pa*abar )...
	/(mu*Pa*Aa) )^(1/(mu-1));

Na_supplied = min(Na_supplied,1);

wR  = mu*Pa*Aa*Na_supplied^(mu-1);

uss	= lambda*(1-Na_supplied)*(1-pQ)/(pQ + lambda*(1-pQ) );

a_u = (be*wc - cbar + alpha/(1-alpha)*Pa*abar )/(Pa)*(1-alpha);
a_e = (wc - cbar + alpha/(1-alpha)*Pa*abar )/(Pa)*(1-alpha);
a_R = (wR - cbar + alpha/(1-alpha)*Pa*abar )/(Pa)*(1-alpha);


ag_demand = uss*a_u + Na_supplied*a_R + (1-uss-Na_supplied)*a_e;

excess(2) = ag_demand - Aa*Na_supplied^mu;
if(Na_supplied>=1) excess(2) = excess(2) - Pa^2; end
	

%excess(3) = be*uss - wcPa(1)*wcPa(3)*(1-uss-Na_supplied);

theeconomy= [Na_supplied,uss,Q,J,VeVu(1),VeVu(2)];
