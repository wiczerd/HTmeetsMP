function [excess,theeconomy] = sol_wcPa_fwd(wcPa,theeconomy_tp1,uNa_tm1)
% computes a point in the economy assuming that the quantities are known in
% the t+1 period 

% if only one output is requested, it gives excess (as in excess demand)
% if more are requested it outputs:
% theeconomy{1} : N_a
% theeconomy{2} : u
% theeconomy{3} : Q
% theeconomy{4} : J
% theeconomy{5} : Ve
% theeconomy{6} : Vu


global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be

wc = wcPa(1);
Pa = wcPa(2);

Jp	= theeconomy_tp1(4);
Vep	= theeconomy_tp1(5);
Vup	= theeconomy_tp1(6);

u_tm1	= uNa_tm1(2);
Na_tm1	= uNa_tm1(1);


J	= (1-lambda)*(Ym-wc + beta*Jp);
Q	= (kappa/(Ym-wc + beta*J)/Amf)^(-1/eta);
pQ	= Amf*Q^(1-eta);
% solve for Ve Vu as a linear system
flowVe	= alpha^alpha*(Pa/(1-alpha)).^(alpha-1).*(wc-cbar-Pa*abar);
flowVu	= alpha^alpha*(Pa/(1-alpha)).^(alpha-1).*(be-cbar-Pa*abar);

VeVucoe	= [0, lambda;...
	0, 0];
% constant values:
VeVucon	= [(1-lambda)*flowVe;pQ*flowVe + (1-pQ)*flowVu];
% continuation values:
VeVucont =  beta*[(1-lambda)*Vep;pQ*Vep + (1-pQ)*Vup];

VeVu	= (eye(2)- VeVucoe)\(VeVucon + VeVucont);

wc_implied = ( alpha^alpha*(Pa/(1-alpha)).^(alpha-1)*...
	( (1-theta)*be + (theta-1)*beta*(Vep - Vup) ) ...
	+ theta*(Ym + beta*Jp) )...
	/(theta+alpha^alpha*(Pa/(1-alpha)).^(alpha-1)*(1-theta)); 

excess(1) = wc_implied - wc;

%solve for the ag economy stuff

Na_supplied = ( ( (VeVu(2)-beta*Vup)/(alpha^alpha*(Pa/(1-alpha))^(alpha-1)) + cbar +Pa*abar )...
	/(mu*Pa*Aa) )^(1/(mu-1));

Na_supplied = min(Na_supplied,1);

wR	= mu*Pa*Aa*Na_supplied^(mu-1);

ut	= (1-pQ)*(u_tm1 + lambda*(1-u_tm1 -Na_supplied) + (Na_tm1 - Na_supplied));

a_u	= (be - cbar + alpha/(1-alpha)*Pa*abar )/(Pa)*(1-alpha);
a_e	= (wc - cbar + alpha/(1-alpha)*Pa*abar )/(Pa)*(1-alpha);
a_R	= (wR - cbar + alpha/(1-alpha)*Pa*abar )/(Pa)*(1-alpha);


ag_demand  = ut*a_u + Na_supplied*a_R + (1-ut-Na_supplied)*a_e;
excess(2)  = ag_demand - Aa*Na_supplied^mu;
if(Na_supplied>=1) excess(2) = excess(2) - Pa^2; end

theeconomy = [Na_supplied,ut,Q,J,VeVu(1),VeVu(2)];
