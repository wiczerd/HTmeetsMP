% This program computes the transition for Harris-Todaro Meets
% Mortensen-Pissarides

cd ~/Documents/CurrResearch/Devt/Computation

global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau

cbar	= -0.6;%-0.6;
abar	= 0.09; %0.012;
Aa	= 1.0;
beta	= 0.99;
eta	= 0.72;
Ym	= 1.12;
lambda	= 0.03;
kappa	= 0.19;
theta	= 0.72;
Amf	= 0.25;
mu	= 0.99;
alpha	= 0.992;
be	= 0.3;
tau	= 0.0;

% util = @(c,a) (c-cbar).^alpha.*(a-abar).^(1-alpha);
% 
% utilc = @(c,a) alpha*((a-abar)./(c-cbar)).^(1-alpha);
% utila = @(c,a) (1-alpha)*((c-cbar)./(a-abar)).^alpha;
% 
% ucua  = @(c,a) alpha/(1-alpha)*((a-abar)./(c-cbar));
% 
% indiru = @(i,Pa) alpha^alpha*(Pa/(1-alpha)).^(alpha-1).*(i-cbar-Pa*abar) ;
% 
% qrt = @(Q) Q.^-eta;
% prt = @(Q) Q.(1-eta);

%% solve for the steady state

sol_wcPa_ss([.5,.4])
%pos_solwcPa = @(logwcPa) sol_wcPa_ss(exp(logwcPa));
pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);


[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)]);
% loop on tau
tauH = .1;tauL=0.;
for itertau = 1:100
	tau = 0.5*tauH+0.5*tauL;
	[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,logssp,optimset('Display','off'));
	wcPa_ss = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
	% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
	[excess_ss,sseconomy] = sol_wcPa_ss([wcPa_ss]);
	budget_def = be*sseconomy(2) - wcPa_ss(1)*tau*(1-sseconomy(2)-sseconomy(1));
	if(abs(budget_def)<1e-6 || (tauH-tauL)<1e-6)
		break;
	elseif (budget_def < 0)
		tauH=tau;
	elseif(budget_def > 0)
		tauL=tau;
	end
end
	
%% calibrate it

% first calibrate abar and Amf so that matches 5% in ag and 5%
% unemployment
[x,fval_cal,exitflag,out] = fminsearch(@cal_ss,log([abar,Amf]));
[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)]);
wcPa_ss = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
[excess_devd,devd_economy] = sol_wcPa_ss(wcPa_ss);

% now calibrate it so that I get 90% in agriculture by chaning Aa
AaH = Aa;  AaL = .01;

for Aaiter = 1:100
	Aa = .5*(AaH + AaL);
	if Aaiter<10
		p0 = [tan(.5*pi/Ym-pi/2) log(.5)];
	else
		p0 = logssp;
	end
	tauH = 0.05; tauL=0.001;
	for itertau = 1:100
		tau = 0.5*tauH+0.5*tauL;
		[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,p0,optimset('Display','off'));
		wcPa_ss = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
		% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
		[excess_undevd,undevd_economy] = sol_wcPa_ss(wcPa_ss);
		budget_def = be*sseconomy(2) - wcPa_ss(1)*tau*(1-sseconomy(2)-sseconomy(1));
		if(abs(budget_def)<1e-6 || (tauH-tauL)<1e-6)
			break;
		elseif (budget_def < 0)
			tauH=tau;
		elseif(budget_def > 0)
			tauL=tau;
		end
	end
	if( abs(undevd_economy(1)-.9) <1e-6 || (AaH-AaL)<1e-6 )
		break;
	elseif(undevd_economy(1)-.9 < 0)
		AaH = Aa;
	elseif(undevd_economy(1)-.9 >0)
		AaL = Aa;
	end
end

%%  Calibrate by moving Amf

Aa = 1.0; Amf = .03;

p0 = [tan(.5*pi/Ym-pi/2) log(.5)];
[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,p0,optimset('Display','off'));
wcPa_ss = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
[lowfnd_excess,lowfnd_economy] = sol_wcPa_ss(wcPa_ss);

%% back one period
%Aa = Aa*.99;
%pos_solwcPa = @(logwcPa) sol_wcPa(exp(logwcPa),sseconomy);
%[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,log(wcPa_ss));
%wcPa_m1 = exp(logssp);
%[excess_m1,economy_m1] = sol_wcPa(wcPa_m1,sseconomy);
