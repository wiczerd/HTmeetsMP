% solves transition paths from low to high steady states

% This program computes the transition for Harris-Todaro Meets
% Mortensen-Pissarides

% This program experiments with slower transitions.

cd ~/Documents/CurrResearch/Devt/Computation

global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau converged

% load everything else (calibration and such):
load trans_AaAmf/calAa_linYmAmf/trans_space_calAa_linYmAmf.mat
%load trans_AaAmf/trans_space_USA.mat 
%solution-level parameters
TT_extra = 200;
TX = TT_extra;
T0 = TT;
TT = T0 + TT_extra;
max_transiter = 41;
save_plots =1;
param_update = 0.5;

C2_prodpaths = 0;
cub_prodpaths = 0;

converged = 0;

%% transition backwards 
%  First, transition using the approximation that unemployment is the steady state.  
%  This doesn't matter explicitly for any of the choices, just market clearing.

Aa = Aa_devd;

%%  set up new sequences

%set up the stunted growth paths in a range, and then put them together for
%the x-section

Ncountry = 10;

Ymfac = ones(Ncountry,1); %change this to be a range
Aafac = ones(Ncountry,1); %change this to be a range
Amffac = ones(Ncountry,1); %change this to be a range

Aafac = linspace(1/Ncountry,Ncountry,Ncountry);

for ci =1:Ncountry

	Ymfac_ci  = Ymfac(ci);
	Aafac_ci  = Aafac(ci);
	Amffac_ci = Amffac(ci);
	if(C2_prodpaths ==1)

		AaTT = (Aa_path(end) - Aa_path(1))/T0/10*T0 + Aa_path(1);
		rt_chng_pwr	= log( (AaTT - Aa_path(1))/(Aa_path(end) - Aa_path(1)) )/log(T0/(TT));
		Aa_path_chng	= linspace(0,1,T0+TT_extra).^rt_chng_pwr;
		Aa_stunted_path		= (1-Aa_path_chng)*Aa_path(1)+Aa_path_chng*Aa_path(end);

		YmTT = (Ym_path(end) - Ym_path(1))/T0/10*T0 + Ym_path(1);
		rt_chng_pwr	= log( (YmTT - Ym_path(1))/(Ym_path(end) - Ym_path(1)) )/log(T0/(T0+TT_extra));
		Ym_path_chng	= linspace(0,1,T0+TT_extra).^rt_chng_pwr;
		Ym_stunted_path		= (1-Ym_path_chng)*Ym_path(1)+Ym_path_chng*Ym_path(end);

		AmfTT = (Amf_path(end) - Amf_path(1))/T0/10*T0 + Amf_path(1);
		rt_chng_pwr	= log( (AmfTT - Amf_path(1))/(Amf_path(end) - Amf_path(1)) )/log(T0/(T0+TT_extra));
		Amf_path_chng	= linspace(0,1,T0+TT_extra).^rt_chng_pwr;
		Amf_stunted_path		= (1-Amf_path_chng)*Amf_path(1)+Amf_path_chng*Amf_path(end);
	elseif(cub_prodpaths == 1)
		MM = [1 0 0 0;1 T0 T0^2 T0^3;1 (T0+TX) (T0+TX)^2 (T0+TX)^3;0 1 2*T0 3*T0^2]; 

		AaTT = (Aa_path(end) - Aa_path(1))/T0/10*T0 + Aa_path(1);
		Mb = [Aa_path(1) AaTT Aa_path(T0) 0]';
		Aa_coef = MM\Mb;
		Aa_stunted_path = [ones(1,T0+TX);1:(T0+TX);(1:(T0+TX)).^2;(1:(T0+TX)).^3]'*Aa_coef;

		YmTT = (Ym_path(end) - Ym_path(1))/T0/10*T0 + Ym_path(1);
		Mb = [Ym_path(1) YmTT Ym_path(T0) 0]';
		Ym_coef = MM\Mb;
		Ym_stunted_path = [ones(1,T0+TX);1:(T0+TX);(1:(T0+TX)).^2;(1:(T0+TX)).^3]'*Ym_coef;

		AmfTT = (Amf_path(end) - Amf_path(1))/T0/10*T0 + Amf_path(1);
		Mb = [Amf_path(1) YmTT Amf_path(T0) 0]';
		Amf_coef = MM\Mb;
		Amf_stunted_path = [ones(1,T0+TX);1:(T0+TX);(1:(T0+TX)).^2;(1:(T0+TX)).^3]'*Amf_coef;


	else
		tmp_growth = log(Aa_path(2:T0) ./ Aa_path(1:T0-1))*Aafac_ci;
		Aa_stunted_path = [Aa_path(1) exp( log(Aa_path(1)) + cumsum(tmp_growth) )];
		% augment it with 200 more periods to get to the final SS
		Aa_stunted_path = [Aa_stunted_path linspace(Aa_stunted_path(T0),Aa_path(T0),TT_extra)];

		tmp_growth = log(Amf_path(2:T0) ./ Amf_path(1:T0-1))*Amffac_ci;
		Amf_stunted_path = [Amf_path(1) exp( log(Amf_path(1)) + cumsum(tmp_growth) )];
		Amf_stunted_path = [Amf_stunted_path linspace(Amf_stunted_path(end),Amf_path(T0),TT_extra)];

		tmp_growth = log(Ym_path(2:T0) ./ Ym_path(1:T0-1))*Ymfac_ci;
		Ym_stunted_path = [Ym_path(1) exp( log(Ym_path(1)) + cumsum(tmp_growth) )];
		Ym_stunted_path = [Ym_stunted_path linspace(Ym_stunted_path(end),Ym_path(T0),TT_extra)];

	end

	Amf_path_ci = [Amf_path ones(1,TT_extra)*Amf_path(T0)];
	Aa_path_ci = [Aa_path  ones(1,TT_extra)*Aa_path(T0)];
	Ym_path_ci = [Ym_path ones(1,TT_extra)*Ym_path(T0)];
	%Amf_stunted_path = Amf_stunted_path(1:TT);
	%Ym_stunted_path = Ym_stunted_path(1:TT);
	%Aa_stunted_path = Aa_stunted_path(1:TT);

	%make this alternate
	%Amf_path_ci=	Amf_stunted_path;
	Aa_path_ci =	Aa_stunted_path;
	%Ym_path_ci =	Ym_stunted_path;
	
	% initially put it to p0 or old price path?
	%price_path  = [p0_trans p0_trans(:,TT)*ones(TT_extra)];
	price_path_ci  = [price_path; ones(TT_extra,1)*price_path(T0,:)];


	% set up storage for things computed over the path
	trans_path  = zeros(TT,size(devd_economy,2));
	excess_path = zeros(TT,2);
	uss_path = zeros(TT,1);
	bad_periods = zeros(TT,1);
	tau_path = zeros(TT,1);

	trans_path(:,2) = -1; % this is an initialization that tells it to replace ut with uss.


	trans_economy	= devd_economy;
	trans_path(TT,:)= trans_economy;

	trans_economy	= undevd_economy;
	trans_path(1,:)= trans_economy;

	trans_path_back = trans_path;
	price_path_fwd = zeros(size(price_path_ci));
	price_path_back= zeros(size(price_path_ci));
	solpath_back = zeros(size(price_path_ci));
	solpath_fwd = zeros(size(price_path_ci));

	Aa_implied_back = ones(TT,1)*Aa_devd;
	Aa_implied_fwd  = ones(TT,1)*Aa_devd;


	price_path_back(TT,:) = wcPa_devd;
	price_path_fwd(TT,:) = wcPa_devd;

	%%

	% run through the transition at least once, just to be sure it converged.
	for trans_iter =1:max_transiter
		logwP = devd_logssp;
		wcPa_t= wcPa_devd;
		for t = TT-1:-1:1
			Aa = Aa_path_ci(t);
			Amf = Amf_path_ci(t);
			Ym = Ym_path_ci(t);


			logwP(2) = log(0.5*price_path_ci(t,2) + 0.5*exp(logwP(2)));
			logwP(1) = tan( (0.5*price_path_ci(t,1)+ 0.5*wcPa_t(1) )*pi/Ym-pi/2);

			pos_solwcPa = @(wcPa) sol_wcPa([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))],trans_path(t+1,:),trans_path(t,2) );
			tauH = 0.05; tauL=0.001;
			for itertau = 1:100
				tau = 0.5*tauH+0.5*tauL;
				p0 = logwP;
				[logwP, fval,exitflag,output,J] = fsolve(pos_solwcPa,p0,optimset('Display','off'));
				if exitflag<1 %try with a new starting point
					p0 = solpath_back(t+1,:);
					[logwP, fval,exitflag,output,J] = fsolve(pos_solwcPa,p0,optimset('Display','off'));
				end
				wcPa_t = [(atan(logwP(1))+pi/2)*Ym/pi exp(logwP(2))];
				% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
				[excess_trans,trans_economy] = sol_wcPa(wcPa_t,trans_path(t+1,:),trans_path(t,2));

				%this has budget clearing in proprotional replacement
				budget_def = be*trans_economy(2) - tau*(1-trans_economy(2)-trans_economy(1));
				%this has budget clearing in absolute replacement
				%budget_def = be*trans_economy(2) - wcAa_t(1)*tau*(1-trans_economy(2)-trans_economy(1));
				if(abs(budget_def)<1e-6 || (tauH-tauL)<1e-6)
					break;
				elseif (budget_def < 0)
					tauH=tau;
				elseif(budget_def > 0)
					tauL=tau;
				end
			end
			tau_path(t) = tau;
			if trans_economy(1) >= 1. || excess_trans*excess_trans' >1e-4;
				trans_economy = trans_path_back(t+1,:);
				bad_periods(t) = 1;
				solpath_back(t,:) = solpath_back(t+1,:);
			else
				% this is the guess for the path forward
				solpath_back(t,:) = logwP;
			end
			trans_path_back(t,:) = trans_economy;
			trans_path(t,:) = trans_economy;
			excess_path(t,:)= excess_trans;
			solpath_back(t,:) = logwP;
			price_path_back(t,:) = [wcPa_t(1) wcPa_t(2)]; 
			pQ	= Amf*trans_economy(3)^(1-eta);
			uss_path(t) = lambda*(1-trans_economy(1))*(1-pQ)/(pQ + lambda*(1-pQ) );
		end

		price_path_ci(:,1) = price_path_back(:,1); % replace wages, but hold fixed the Pa;
		price_path_fwd(1,:) = wcPa_t;

		% now go forward to get quantities right
		for t = 2:TT-1
			Aa = Aa_path_ci(t);
			Amf = Amf_path_ci(t);
			Ym = Ym_path_ci(t);

			pos_solwcPa = @(wcPa) sol_wcPa_fwd([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))],trans_path(t+1,:),trans_path(t-1,1:2));

			tau = tau_path(t);

			p0 = solpath_back(t,:); 
			[logwP, fval,exitflag,output,J] = fsolve(pos_solwcPa,p0,optimset('Display','off'));
			wcPa_t = [(atan(logwP(1))+pi/2)*Ym/pi exp(logwP(2))];
			if(exitflag<=0)
				p0 = solpath_fwd(t-1,:);
				[logwP, fval,exitflag,output,J] = fsolve(pos_solwcPa,p0,optimset('Display','off'));
				wcPa_t = [(atan(logwP(1))+pi/2)*Ym/pi exp(logwP(2))];
			end
			if(exitflag<=0)
				logwP = p0;
				wcPa_t = [price_path_ci(t,1)  price_path_ci(t,2)];
			end
			solpath_fwd(t,:) = logwP;

			% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
			[excess_trans,trans_economy] = sol_wcPa_fwd(wcPa_t,trans_path(t+1,:),trans_path(t-1,1:2));

			budget_def = be*trans_economy(2) - tau*(1-trans_economy(2)-trans_economy(1));

			if trans_economy(1) >= 1. || exitflag<0;
				trans_economy = trans_path(t-1,:);
				bad_periods(t) = 1;
			end


			trans_path(t,:) = trans_economy;
			excess_path(t,:)= excess_trans;
			price_path_fwd(t,:) = [wcPa_t(1) wcPa_t(2)]; % or if need to go slower: wcPa_t*.25 + .75*price_path_back(t,:);
		end

		%post transition checks

		trans_path(TT,:) = devd_economy;
		price_dif = abs(price_path_fwd - price_path_back);
		price_path_ci = price_path_fwd;

		resid = sum(price_dif);

		if(resid < 1e-4)
			break;
		end

	end


	%% Compute paths for revenue per worker in Ag & Urban at P_t and P_0

	% rev_pt = [price_path(:,2).*Aa_path'.*trans_path(:,1).^mu Ym_path'];  % ag,man
	% percaprev_pt = [price_path(:,2).*Aa_path'.*trans_path(:,1).^(mu-1) Ym_path'./(1-trans_path(:,1)-trans_path(:,2) )];  % ag,man
	% rev_p0 = [price_path(1,2).*Aa_path'.*trans_path(:,1).^mu Ym_path'];  % ag,man
	% percaprev_p0 = [price_path(1,2).*Aa_path'.*trans_path(:,1).^(mu-1) Ym_path'./(1-trans_path(:,1)-trans_path(:,2) )];  % ag,man
	% % this is using the TX period prices: not fully developed
	% rev_pTX = [price_path(TT-TX,2).*Aa_path'.*trans_path(:,1).^mu Ym_path'];  % ag,man
	% percaprev_pTX = [price_path(TT-TX,2).*Aa_path'.*trans_path(:,1).^(mu-1) Ym_path'./(1-trans_path(:,1)-trans_path(:,2) )];  % ag,man
	% % this is using the TT period prices: fully developed
	% rev_pTT = [price_path(TT,2).*Aa_path'.*trans_path(:,1).^mu Ym_path'];  % ag,man
	% percaprev_pTT = [price_path(TT,2).*Aa_path'.*trans_path(:,1).^(mu-1) Ym_path'./(1-trans_path(:,1)-trans_path(:,2) )];  % ag,man
	% 
	rev_pt = [price_path_ci(:,2).*Aa_path_ci'.*trans_path(:,1).^mu Ym_path_ci'.*(1-trans_path(:,2)-trans_path(:,1))];  % ag,man
	percaprev_pt = [price_path_ci(:,2).*Aa_path_ci'.*trans_path(:,1).^(mu-1) Ym_path_ci'./(1-trans_path(:,2) )];  % ag,man
	rev_p0 = [price_path_ci(1,2).*Aa_path_ci'.*trans_path(:,1).^mu Ym_path_ci'.*(1-trans_path(:,2)-trans_path(:,1))];  % ag,man
	percaprev_p0 = [price_path_ci(1,2).*Aa_path_ci'.*trans_path(:,1).^(mu-1) Ym_path_ci'./(1-trans_path(:,2) )];  % ag,man
	% this is using the TX period prices: not fully developed
	rev_pTX = [price_path_ci(TT-TX,2).*Aa_path_ci'.*trans_path(:,1).^mu Ym_path_ci'.*(1-trans_path(:,2)-trans_path(:,1))];  % ag,man
	percaprev_pTX = [price_path_ci(TT-TX,2).*Aa_path_ci'.*trans_path(:,1).^(mu-1) Ym_path_ci'./(1-trans_path(:,2) )];  % ag,man
	% this is using the TT period prices: fully developed
	rev_pTT = [price_path_ci(TT,2).*Aa_path_ci'.*trans_path(:,1).^mu Ym_path_ci'.*(1-trans_path(:,2)-trans_path(:,1))];  % ag,man
	percaprev_pTT = [price_path_ci(TT,2).*Aa_path_ci'.*trans_path(:,1).^(mu-1) Ym_path_ci'./(1-trans_path(:,2) )];  % ag,man



	upath = trans_path(:,2)./(1-trans_path(:,1));
	mpath = (-trans_path(2:TT,1) + trans_path(1:TT-1,1))./trans_path(1:TT-1,1);

	%%
	cd trans_AaAmf/xsec_AaAmfYm
	save(['trans_space_' num2str(ci) '.mat']);

	cd ~/Documents/CurrResearch/Devt/Computation


end % loop over country

apg_ci = zeros(TT,Ncountry);
u_ci = zeros(TT,Ncountry);
y_ci = zeros(TT,Ncountry);

for ci =1:Ncountry
	load(['trans_AaAmf/xsec_AaAmfYm/trans_space_' num2str(ci) '.mat']);
	apg_ci(:,ci)  = percaprev_pt(:,2)./percaprev_pt(:,1);
	u_ci(:,ci)  = upath;
	y_ci(:,ci)  = log( trans_path(:,1).*(Aa_path_ci'.*trans_path(:,1).^alpha)+(1. - trans_path(:,1) - trans_path(:,2)).*Ym_path_ci' );
end


TO = TT- TT_extra;
tt = 10;

cd trans_AaAmf/xsec_AaAmfYm;

figure(1);
scatter(reshape(y_ci([1:TO],:),T0*Ncountry,1),reshape(u_ci([1:T0],:),T0*Ncountry,1),'LineWidth',2);
set(gcf,'color','white');
title('Unemployment Scatter - All Perdiods');
grid on;
if (save_plots == 1) saveas(gca,'u_xsec','eps2c'); end
if (save_plots == 1) saveas(gca,'u_xsec.png'); end

figure(2);
scatter(y_ci(1,:),u_ci(1,:),'LineWidth',2);
set(gcf,'color','white');
title('Unemployment Scatter - Perdiod 1');
grid on;
if (save_plots == 1) saveas(gca,'u_xsec1','eps2c'); end
if (save_plots == 1) saveas(gca,'u_xsec1.png'); end

figure(3);
scatter(y_ci(T0,:),u_ci(T0,:),'LineWidth',2);
set(gcf,'color','white');
title('Unemployment Scatter - Perdiod T');
grid on;
if (save_plots == 1) saveas(gca,'u_xsecT','eps2c'); end
if (save_plots == 1) saveas(gca,'u_xsecT.png'); end

figure(4);
scatter(y_ci(tt,:),u_ci(tt,:),'LineWidth',2);
set(gcf,'color','white');
title(['Unemployment Scatter - Perdiod t=' num2str(tt)]);
grid on;
if (save_plots == 1) saveas(gca,'u_xsect','eps2c'); end
if (save_plots == 1) saveas(gca,'u_xsect.png'); end

figure(11);
scatter(reshape(y_ci([1:TO],:),T0*Ncountry,1),reshape(apg_ci([1:T0],:),T0*Ncountry,1),'LineWidth',2);
set(gcf,'color','white');
title('APG Scatter - All Perdiods');
grid on;
if (save_plots == 1) saveas(gca,'apg_xsec','eps2c'); end
if (save_plots == 1) saveas(gca,'apg_xsec.png'); end

figure(12);
scatter(y_ci(1,:),apg_ci(1,:),'LineWidth',2);
set(gcf,'color','white');
title('APG Scatter - Perdiod 1');
grid on;
if (save_plots == 1) saveas(gca,'apg_xsec1','eps2c'); end
if (save_plots == 1) saveas(gca,'apg_xsec1.png'); end

figure(13);
scatter(y_ci(T0,:),apg_ci(T0,:),'LineWidth',2);
set(gcf,'color','white');
title('APG Scatter - Perdiod T');
grid on;
if (save_plots == 1) saveas(gca,'apg_xsecT','eps2c'); end
if (save_plots == 1) saveas(gca,'apg_xsecT.png'); end

figure(14);
scatter(y_ci(tt,:),apg_ci(tt,:),'LineWidth',2);
set(gcf,'color','white');
title(['APG Scatter - Perdiod t=' num2str(tt)]);
grid on;
if (save_plots == 1) saveas(gca,'apg_xsect','eps2c'); end
if (save_plots == 1) saveas(gca,'apg_xsect.png'); end