% solves transition paths from low to high steady states

% This program computes the transition for Harris-Todaro Meets
% Mortensen-Pissarides

% This program experiments with slower transitions.

cd ~/Documents/CurrResearch/Devt/Computation

global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau converged

% load everything else (calibration and such):
load trans_AaAmf/calAa_linYmAmf/trans_space_calAa_linYmAmf.mat

%solution-level parameters
TT_extra = 200;
TX = TT_extra;
max_transiter = 41;
save_plots =1;
param_update = 0.5;
Amf_stunted= 0;
Ym_stunted = 1;
Aa_stunted = 0;

C2_prodpaths = 0;
cub_prodpaths = 0;

converged = 0;

%% transition backwards 
%  First, transition using the approximation that unemployment is the steady state.  
%  This doesn't matter explicitly for any of the choices, just market clearing.

Aa = Aa_devd;

%%  set up new sequences

%set up the stunted growth paths

if(C2_prodpaths ==1)

	AaTT = (Aa_path(end) - Aa_path(1))/TT/10*TT + Aa_path(1);
	rt_chng_pwr	= log( (AaTT - Aa_path(1))/(Aa_path(end) - Aa_path(1)) )/log(TT/(TT+TT_extra));
	Aa_path_chng	= linspace(0,1,TT+TT_extra).^rt_chng_pwr;
	Aa_stunted_path		= (1-Aa_path_chng)*Aa_path(1)+Aa_path_chng*Aa_path(end);

	YmTT = (Ym_path(end) - Ym_path(1))/TT/10*TT + Ym_path(1);
	rt_chng_pwr	= log( (YmTT - Ym_path(1))/(Ym_path(end) - Ym_path(1)) )/log(TT/(TT+TT_extra));
	Ym_path_chng	= linspace(0,1,TT+TT_extra).^rt_chng_pwr;
	Ym_stunted_path		= (1-Ym_path_chng)*Ym_path(1)+Ym_path_chng*Ym_path(end);

	AmfTT = (Amf_path(end) - Amf_path(1))/TT/10*TT + Amf_path(1);
	rt_chng_pwr	= log( (AmfTT - Amf_path(1))/(Amf_path(end) - Amf_path(1)) )/log(TT/(TT+TT_extra));
	Amf_path_chng	= linspace(0,1,TT+TT_extra).^rt_chng_pwr;
	Amf_stunted_path		= (1-Amf_path_chng)*Amf_path(1)+Amf_path_chng*Amf_path(end);
elseif(cub_prodpaths == 1)
	MM = [1 0 0 0;1 TT TT^2 TT^3;1 (TT+TX) (TT+TX)^2 (TT+TX)^3;0 1 2*TT 3*TT^2]; 
	
	AaTT = (Aa_path(end) - Aa_path(1))/TT/10*TT + Aa_path(1);
	Mb = [Aa_path(1) AaTT Aa_path(TT) 0]';
	Aa_coef = MM\Mb;
	Aa_stunted_path = [ones(1,TT+TX);1:(TT+TX);(1:(TT+TX)).^2;(1:(TT+TX)).^3]'*Aa_coef;
	
	YmTT = (Ym_path(end) - Ym_path(1))/TT/10*TT + Ym_path(1);
	Mb = [Ym_path(1) YmTT Ym_path(TT) 0]';
	Ym_coef = MM\Mb;
	Ym_stunted_path = [ones(1,TT+TX);1:(TT+TX);(1:(TT+TX)).^2;(1:(TT+TX)).^3]'*Ym_coef;
	
	AmfTT = (Amf_path(end) - Amf_path(1))/TT/10*TT + Amf_path(1);
	Mb = [Amf_path(1) YmTT Amf_path(TT) 0]';
	Amf_coef = MM\Mb;
	Amf_stunted_path = [ones(1,TT+TX);1:(TT+TX);(1:(TT+TX)).^2;(1:(TT+TX)).^3]'*Amf_coef;
	
	
else
	tmp_growth = log(Aa_path(2:end) ./ Aa_path(1:end-1))/10;
	Aa_stunted_path = [Aa_path(1) exp( log(Aa_path(1)) + cumsum(tmp_growth) )];
	% augment it with 200 more periods to get to the final SS
	Aa_stunted_path = [Aa_stunted_path linspace(Aa_stunted_path(end),Aa_path(end),TT_extra)];

	tmp_growth = log(Amf_path(2:end) ./ Amf_path(1:end-1))/10;
	Amf_stunted_path = [Amf_path(1) exp( log(Amf_path(1)) + cumsum(tmp_growth) )];
	Amf_stunted_path = [Amf_stunted_path linspace(Amf_stunted_path(end),Amf_path(end),TT_extra)];

	tmp_growth = log(Ym_path(2:end) ./ Ym_path(1:end-1))/10;
	Ym_stunted_path = [Ym_path(1) exp( log(Ym_path(1)) + cumsum(tmp_growth) )];
	Ym_stunted_path = [Ym_stunted_path linspace(Ym_stunted_path(end),Ym_path(end),TT_extra)];

end

Amf_path = [Amf_path ones(1,TT_extra)*Amf_path(TT)];
Aa_path = [Aa_path  ones(1,TT_extra)*Aa_path(TT)];
Ym_path = [Ym_path ones(1,TT_extra)*Ym_path(TT)];
%Amf_stunted_path = Amf_stunted_path(1:TT);
%Ym_stunted_path = Ym_stunted_path(1:TT);
%Aa_stunted_path = Aa_stunted_path(1:TT);

if Amf_stunted==1;	Amf_path=	Amf_stunted_path;	end
if Aa_stunted== 1;	Aa_path =	Aa_stunted_path;	end
if Ym_stunted== 1;	Ym_path =	Ym_stunted_path;	end


% initially put it to p0 or old price path?
%price_path  = [p0_trans p0_trans(:,TT)*ones(TT_extra)];
price_path  = [price_path; ones(TT_extra,1)*price_path(TT,:)];

% allocate all of the tracking matrices so that it has the new size
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TT = TT + TT_extra;

% set up storage for things computed over the path
trans_path  = zeros(TT,size(devd_economy,2));
excess_path = zeros(TT,2);
uss_path = zeros(TT,1);

trans_path(:,2) = -1; % this is an initialization that tells it to replace ut with uss.


trans_economy	= devd_economy;
trans_path(TT,:)= trans_economy;

trans_economy	= undevd_economy;
trans_path(1,:)= trans_economy;

trans_path_back = trans_path;
price_path_fwd = zeros(size(price_path));
price_path_back= zeros(size(price_path));
solpath_back = zeros(size(price_path));

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
		Aa = Aa_path(t);
		Amf = Amf_path(t);
		Ym = Ym_path(t);

		
		logwP(2) = log(0.5*price_path(t,2) + 0.5*exp(logwP(2)));
		logwP(1) = tan( (0.5*price_path(t,1)+ 0.5*wcPa_t(1) )*pi/Ym-pi/2);

		pos_solwcPa = @(wcPa) sol_wcPa([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))],trans_path(t+1,:),trans_path(t,2) );
		tauH = 0.05; tauL=0.001;
		for itertau = 1:100
			tau = 0.5*tauH+0.5*tauL;
			p0 = logwP;
			[logwP, fval,exitflag,output,J] = fsolve(pos_solwcPa,p0,optimset('Display','off'));
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
		trans_path_back(t,:) = trans_economy;
		trans_path(t,:) = trans_economy;
		excess_path(t,:)= excess_trans;
		solpath_back(t,:) = logwP;
		price_path_back(t,:) = [wcPa_t(1) wcPa_t(2)]; 
		pQ	= Amf*trans_economy(3)^(1-eta);
		uss_path(t) = lambda*(1-trans_economy(1))*(1-pQ)/(pQ + lambda*(1-pQ) );
	end

	price_path(:,1) = price_path_back(:,1); % replace wages, but hold fixed the Pa;
	price_path_fwd(1,:) = wcPa_t;

	% now go forward to get quantities right
	for t = 2:TT-1
		Aa = Aa_path(t);
		Amf = Amf_path(t);
		Ym = Ym_path(t);
				
		pos_solwcPa = @(wcPa) sol_wcPa_fwd([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))],trans_path(t+1,:),trans_path(t-1,1:2));

		tauH = 0.05; tauL=0.001;
		for itertau = 1:100
			tau = 0.5*tauH+0.5*tauL;

			p0 = solpath_back(t,:); 
			[logwP, fval,exitflag,output,J] = fsolve(pos_solwcPa,p0,optimset('Display','off'));
			wcPa_t = [(atan(logwP(1))+pi/2)*Ym/pi exp(logwP(2))];
			% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
			if(exitflag<=0)
				logwP = p0;
				wcPa_t = [price_path(t,1)  price_path(t,2)];
			end

			[excess_trans,trans_economy] = sol_wcPa_fwd(wcPa_t,trans_path(t+1,:),trans_path(t-1,1:2));

			budget_def = be*trans_economy(2) - tau*(1-trans_economy(2)-trans_economy(1));
			if(abs(budget_def)<1e-6 || (tauH-tauL)<1e-6)
				break;
			elseif (budget_def < 0)
				tauH=tau;
			elseif(budget_def > 0)
				tauL=tau;
			end
		end
		trans_path(t,:) = trans_economy;
		excess_path(t,:)= excess_trans;
		price_path_fwd(t,:) = [wcPa_t(1) wcPa_t(2)]; % or if need to go slower: wcPa_t*.25 + .75*price_path_back(t,:);
	end
	
	%post transition checks
	
	trans_path(TT,:) = devd_economy;
	price_dif = abs(price_path_fwd - price_path_back);
	price_path = price_path_fwd;

	resid = sum(price_dif);

	if(resid < 1e-4)
		break;
	end
	
end


%% Compute paths for revenue per worker in Ag & Urban at P_t and P_0

rev_pt = [price_path(:,2).*Aa_path'.*trans_path(:,1).^mu Ym_path'];  % ag,man
percaprev_pt = [price_path(:,2).*Aa_path'.*trans_path(:,1).^(mu-1) Ym_path'./(1-trans_path(:,1))];  % ag,man
rev_p0 = [price_path(1,2).*Aa_path'.*trans_path(:,1).^mu Ym_path'];  % ag,man
percaprev_p0 = [price_path(1,2).*Aa_path'.*trans_path(:,1).^(mu-1) Ym_path'./(1-trans_path(:,1))];  % ag,man
% this is using the TX period prices: not fully developed
rev_pTX = [price_path(TT-TX,2).*Aa_path'.*trans_path(:,1).^mu Ym_path'];  % ag,man
percaprev_pTX = [price_path(TT-TX,2).*Aa_path'.*trans_path(:,1).^(mu-1) Ym_path'./(1-trans_path(:,1))];  % ag,man
% this is using the TT period prices: fully developed
rev_pTT = [price_path(TT,2).*Aa_path'.*trans_path(:,1).^mu Ym_path'];  % ag,man
percaprev_pTT = [price_path(TT,2).*Aa_path'.*trans_path(:,1).^(mu-1) Ym_path'./(1-trans_path(:,1))];  % ag,man



%%
upath = trans_path(:,2)./(1-trans_path(:,1));
mpath = (-trans_path(2:TT,1) + trans_path(1:TT-1,1))./trans_path(1:TT-1,1);

if (save_plots==1 && Aa_stunted == 1 && Ym_stunted==1) 
	cd trans_AaAmf/calAa_haltYm_haltAa
	save trans_space_calAa_haltYm_haltAa;
elseif (save_plots==1 && Amf_stunted == 1 && Ym_stunted==1) 
	cd trans_AaAmf/calAa_haltYm_haltAmf
	save trans_space_calAa_haltYm_haltAmf;	
elseif (save_plots==1 && Amf_stunted == 1) 
	cd trans_AaAmf/calAa_linYm_haltAmf
	save trans_space_calAa_linYm_haltAmf;
elseif (save_plots==1 && Aa_stunted == 1) 
	cd trans_AaAmf/calAa_linYm_haltAa
	save trans_space_calAa_linYm_haltAa;	
elseif (save_plots==1 && Ym_stunted == 1) 
	cd trans_AaAmf/calAa_linAmf_haltYm
	save trans_space_calAa_linAmf_haltYm;	
end

%
figure(1);
[ax,h1,h2]=plotyy([1:TT],trans_path(:,1),[1:TT],Aa_path);title('Fraction in agriculture');
set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','N_a','A_a');
set(gcf,'color','white');
grid on;
if (save_plots == 1) saveas(gca,'Natrans','eps2c'); end

figure(2);
[ax,h1,h2]=plotyy([2:TT],mpath,[1:TT],Aa_path);title('Rural -> urban rate');
set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','m/N_a','A_a');
set(gcf,'color','white');
grid on;
if (save_plots == 1) saveas(gca,'mtrans','eps2c'); end

figure(3);
[ax,h1,h2]=plotyy([1:TT],upath,[1:TT],Amf_path);title('Unemployment rate');
set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','u','A_m');
set(gcf,'color','white');
grid on;
if (save_plots == 1) saveas(gca,'utrans','eps2c'); end

figure(4);
[ax,h1,h2]=plotyy([1:TT],Amf*trans_path(:,3).^(1-eta),[1:TT],Amf_path);title('Job finding rate');
set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','p(Q)','A_m');
set(gcf,'color','white');
grid on;
if (save_plots == 1) saveas(gca,'pQtrans','eps2c'); end
	
figure(5);
[ax,h1,h2]=plotyy([1:TT],percaprev_pt(:,1)./percaprev_pt(:,2) ,[1:TT], percaprev_pt(:,1));title('Relative Revenue Per Capita, P_t');
set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','NorthEast','P_t y_a/y_m','P_t y_a');
y20 = min(percaprev_pt(:,1));y2Y=max(percaprev_pt(:,1));
y10 = min(percaprev_pt(:,1)./percaprev_pt(:,2));y1Y=max(percaprev_pt(:,1)./percaprev_pt(:,2));
set(ax(1),'YTick',round(([0:6]*(y1Y-y10)/6 +y10)*100)/100 );set(ax(2),'YTick',round(([0:6]*(y2Y-y20)/6 +y20)*100)/100);
set(gcf,'color','white');
grid on;
if (save_plots == 1) saveas(gca,'rev_Pt_trans','eps2c'); end

figure(6);
h=plot([1:TT],percaprev_p0(:,1)./percaprev_p0(:,2) ,[1:TT], percaprev_p0(:,1));title('Relative Revenue Per Capita, P_0');
set(h,'LineWidth',2);legend('Location','North','P_0 y_a/y_m','P_0 y_a');
set(gcf,'color','white');
grid on;
if (save_plots == 1) saveas(gca,'rev_P0_trans','eps2c'); end

figure(7);
h=plot([1:TT],percaprev_pTT(:,1)./percaprev_pTT(:,2) ,[1:TT], percaprev_pTT(:,1));title('Relative Revenue Per Capita, P_T');
set(h,'LineWidth',2);legend('Location','North','P_T y_a/y_m','P_T y_a');
set(gcf,'color','white');
grid on;
if (save_plots == 1) saveas(gca,'rev_PTT_trans','eps2c'); end

figure(8);
[ax,h1,h2]=plotyy([1:TT],price_path(:,2),[1:TT],Aa_path);title('Relative price of agricultural good');
set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','P_t','A_a');
set(gcf,'color','white');
grid on;
if (save_plots == 1) saveas(gca,'agPrice','eps2c'); end


cd ~/Documents/CurrResearch/Devt/Computation


%%
	
load trans_AaAmf/calAa_linYmAmf/trans_space_calAa_linYmAmf;
	upath_calAa_linYmAmf= upath;
	apg_calAa_linYmAmf= percaprev_pt(:,1)./percaprev_pt(:,2);
	apgTT_calAa_linYmAmf= percaprev_pTT(:,1)./percaprev_pTT(:,2);
	gdp_calAa_linYmAmf= log( trans_path(:,1).*(Aa_path'.*trans_path(:,1).^alpha)+(1. - trans_path(:,1)).*Ym_path' );
	relgdp_calAa_linYmAmf= (gdp_calAa_linYmAmf-gdp_calAa_linYmAmf(end));
	Aa_path_calAa_linYmAmf= Aa_path;
	Amf_path_calAa_linYmAmf = Amf_path;
	Ym_path_calAa_linYmAmf = Ym_path;

load trans_AaAmf/calAa_haltYm_haltAa/trans_space_calAa_haltYm_haltAa;
	upath_calAa_haltYm_haltAa = upath;
	apg_calAa_haltYm_haltAa = percaprev_pt(:,1)./percaprev_pt(:,2);
	apgTT_calAa_haltYm_haltAa = percaprev_pTT(:,1)./percaprev_pTT(:,2);
	gdp_calAa_haltYm_haltAa = log( trans_path(:,1).*(Aa_path'.*trans_path(:,1).^alpha)+(1. - trans_path(:,1)).*Ym_path' );
	relgdp_calAa_haltYm_haltAa = (gdp_calAa_haltYm_haltAa-gdp_calAa_linYmAmf(end));
	Aa_path_calAa_haltYm_haltAa  = Aa_path;
	Amf_path_calAa_haltYm_haltAa  = Amf_path;
	Ym_path_calAa_haltYm_haltAa  = Ym_path;

load trans_AaAmf/calAa_linYm_haltAa/trans_space_calAa_linYm_haltAa;
	upath_calAa_linYm_haltAa = upath;
	apg_calAa_linYm_haltAa = percaprev_pt(:,1)./percaprev_pt(:,2);
	apgTT_calAa_linYm_haltAa = percaprev_pTT(:,1)./percaprev_pTT(:,2);
	gdp_calAa_linYm_haltAa = log( trans_path(:,1).*(Aa_path'.*trans_path(:,1).^alpha)+(1. - trans_path(:,1)).*Ym_path' );
	relgdp_calAa_linYm_haltAa = (gdp_calAa_linYm_haltAa - gdp_calAa_linYmAmf(end));
	Aa_path_calAa_linYm_haltAa = Aa_path;
	Amf_path_calAa_linYm_haltAa = Amf_path;
	Ym_path_calAa_linYm_haltAa = Ym_path;

load trans_AaAmf/calAa_linAmf_haltYm/trans_space_calAa_linAmf_haltYm;
	upath_calAa_linAmf_haltYm = upath;
	apg_calAa_linAmf_haltYm = percaprev_pt(:,1)./percaprev_pt(:,2);
	apgTT_calAa_linAmf_haltYm = percaprev_pTT(:,1)./percaprev_pTT(:,2);
	gdp_calAa_linAmf_haltYm = log( trans_path(:,1).*(Aa_path'.*trans_path(:,1).^alpha)+(1. - trans_path(:,1)).*Ym_path' );
	relgdp_calAa_linAmf_haltYm = (gdp_calAa_linAmf_haltYm-gdp_calAa_linYmAmf(end));
	Aa_path_calAa_linAmf_haltYm = Aa_path;
	Amf_path_calAa_linAmf_haltYm = Amf_path;
	Ym_path_calAa_linAmf_haltYm = Ym_path;

save_plots=1;
TO = TT- TT_extra;
figure(10);
plot([1:TO],upath_calAa_linYmAmf(1:TO),'k',[1:TO],upath_calAa_linAmf_haltYm(1:TO),'b', [1:TO], upath_calAa_linYm_haltAa(1:TO),'g',[1:TO], upath_calAa_haltYm_haltAa(1:TO),'r','LineWidth',2);
set(gcf,'color','white');
title('Unemployment paths');
legend('Location','NorthEast','Baseline','Slow y_m','Slow A_a','Slow A_a & y_m');
grid on;
if (save_plots == 1) saveas(gca,'upath_compare','eps2c'); end
if (save_plots == 1) saveas(gca,'upath_compare.png'); end

figure(11);
plot([1:TO],apg_calAa_linYmAmf(1:TO),'k', [1:TO],apg_calAa_linAmf_haltYm(1:TO),'b', [1:TO], apg_calAa_linYm_haltAa(1:TO),'g',[1:TO], apg_calAa_haltYm_haltAa(1:TO),'r','LineWidth',2);
set(gcf,'color','white');
title('Agricultural Productivity Gap, current prices');
legend('Location','NorthWest','Baseline','Slow y_m','Slow A_a','Slow A_a & y_m');
grid on;
if (save_plots == 1) saveas(gca,'apgpath_compare','eps2c'); end
if (save_plots == 1) saveas(gca,'apgpath_compare.png'); end

figure(12);
plot([1:TO],apgTT_calAa_linYmAmf(1:TO),'k', [1:TO],apgTT_calAa_linAmf_haltYm(1:TO),'b', [1:TO], apgTT_calAa_linYm_haltAa(1:TO),'g',[1:TO], apgTT_calAa_haltYm_haltAa(1:TO),'r','LineWidth',2);
set(gcf,'color','white');
title('Agricultural Productivity Gap, Devd prices');
legend('Location','NorthWest','Baseline','Slow y_m','Slow A_a','Slow A_a & y_m');
grid on;
if (save_plots == 1) saveas(gca,'apgDevpath_compare','eps2c'); end
if (save_plots == 1) saveas(gca,'apgDevpath_compare.png'); end


figure(13);
plot([1:TO],Aa_path_calAa_linYmAmf(1:TO),'k:', [1:TO],Aa_path_calAa_linAmf_haltYm(1:TO),'b:', [1:TO], Aa_path_calAa_linYm_haltAa(1:TO),'g:', [1:TO], Aa_path_calAa_haltYm_haltAa(1:TO),'r:', ...
	[1:TO],Ym_path_calAa_linYmAmf(1:TO),'k-', [1:TO],Ym_path_calAa_linAmf_haltYm(1:TO),'b-', [1:TO], Ym_path_calAa_linYm_haltAa(1:TO),'g-', [1:TO], Ym_path_calAa_haltYm_haltAa(1:TO),'r-', ...
	'LineWidth',2);
set(gcf,'color','white');
title('Productivity Paths');
legend('Location','NorthWest','A_a Baseline','A_a - Slow y_m','A_a - Slow A_a','A_a - Slow A_a & y_m','y_m - Baseline','y_m - Slow y_m','y_m - Slow A_a','y_m - Slow A_a & y_m');
grid on;
if (save_plots == 1) saveas(gca,'prodpath_compare','eps2c'); end
if (save_plots == 1) saveas(gca,'prodpath_compare.png'); end


figure(13);
plot(relgdp_calAa_linYmAmf(1:TO),upath_calAa_linYmAmf(1:TO),'k',relgdp_calAa_linAmf_haltYm(1:TO),upath_calAa_linAmf_haltYm(1:TO),'b', relgdp_calAa_linYm_haltAa(1:TO), upath_calAa_linYm_haltAa(1:TO),'g',relgdp_calAa_haltYm_haltAa(1:TO), upath_calAa_haltYm_haltAa(1:TO),'r','LineWidth',2);
set(gcf,'color','white');
title('Unemployment paths');xlabel('log GDP, deviation from devd');
legend('Location','NorthEast','Baseline','Slow y_m','Slow A_a','Slow A_a & y_m');
grid on;
if (save_plots == 1) saveas(gca,'upath_gdp_compare','eps2c'); end
if (save_plots == 1) saveas(gca,'upath_gdp_compare.png'); end

figure(14);
plot(relgdp_calAa_linYmAmf(1:TO),apg_calAa_linYmAmf(1:TO),'k', relgdp_calAa_linAmf_haltYm(1:TO),apg_calAa_linAmf_haltYm(1:TO),'b', relgdp_calAa_linYm_haltAa(1:TO), apg_calAa_linYm_haltAa(1:TO),'g',relgdp_calAa_haltYm_haltAa(1:TO), apg_calAa_haltYm_haltAa(1:TO),'r','LineWidth',2);
set(gcf,'color','white');
title('Agricultural Productivity Gap, current prices');xlabel('log GDP, deviation from devd');
legend('Location','SouthEast','Baseline','Slow y_m','Slow A_a','Slow A_a & y_m');
grid on;
if (save_plots == 1) saveas(gca,'apgpath_gdp_compare','eps2c'); end
if (save_plots == 1) saveas(gca,'apgpath_gdp_compare.png'); end

figure(15);
plot(relgdp_calAa_linYmAmf(1:TO),apgTT_calAa_linYmAmf(1:TO),'k', relgdp_calAa_linAmf_haltYm(1:TO),apgTT_calAa_linAmf_haltYm(1:TO),'b', relgdp_calAa_linYm_haltAa(1:TO), apgTT_calAa_linYm_haltAa(1:TO),'g',relgdp_calAa_haltYm_haltAa(1:TO), apgTT_calAa_haltYm_haltAa(1:TO),'r','LineWidth',2);
set(gcf,'color','white');
title('Agricultural Productivity Gap, Devd prices');xlabel('log GDP, deviation from devd');
legend('Location','NorthEast','Baseline','Slow y_m','Slow A_a','Slow A_a & y_m');
grid on;
if (save_plots == 1) saveas(gca,'apgDevpath_gdp_compare','eps2c'); end
if (save_plots == 1) saveas(gca,'apgDevpath_gdp_compare.png'); end


figure(16);
plot(relgdp_calAa_linYmAmf(1:TO),Aa_path_calAa_linYmAmf(1:TO),'k:', relgdp_calAa_linAmf_haltYm(1:TO),Aa_path_calAa_linAmf_haltYm(1:TO),'b:', relgdp_calAa_linYm_haltAa(1:TO), Aa_path_calAa_linYm_haltAa(1:TO),'g:', relgdp_calAa_haltYm_haltAa(1:TO), Aa_path_calAa_haltYm_haltAa(1:TO),'r:', ...
	 relgdp_calAa_linYmAmf(1:TO),Ym_path_calAa_linYmAmf(1:TO),'k-', relgdp_calAa_linAmf_haltYm(1:TO),Ym_path_calAa_linAmf_haltYm(1:TO),'b-', relgdp_calAa_linYm_haltAa(1:TO), Ym_path_calAa_linYm_haltAa(1:TO),'g-', relgdp_calAa_haltYm_haltAa(1:TO), Ym_path_calAa_haltYm_haltAa(1:TO),'r-', ...
	'LineWidth',2);
set(gcf,'color','white');
title('Productivity Paths');xlabel('log GDP, deviation from devd');
legend('Location','SouthEast','A_a Baseline','A_a - Slow y_m','A_a - Slow A_a','A_a - Slow A_a & y_m','y_m - Baseline','y_m - Slow y_m','y_m - Slow A_a','y_m - Slow A_a & y_m');
grid on;
if (save_plots == 1) saveas(gca,'prodpath_gdp_compare','eps2c'); end
if (save_plots == 1) saveas(gca,'prodpath_gdp_compare.png'); end
