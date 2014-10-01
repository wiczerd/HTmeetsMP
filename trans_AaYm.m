% solves transition paths from low to high steady states

% This program computes the transition for Harris-Todaro Meets
% Mortensen-Pissarides

cd ~/Documents/CurrResearch/Devt/Computation

global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau

TT_data = 150*4;
TT	= 50;%TT_data;
save_plots =0;

cbar	= 0;%-0.6; % note this is the inverse because I changed the util function
abar	= 0.2; % this gets changed below in the calibration
Aa	= 1.0;
beta	= 0.99;
eta	= 0.72;
Ym	= 1.12;
lambda	= 0.03;
kappa	= 0.19;
theta	= 0.72;
Amf	= 0.5;
mu	= 0.99;
alpha	= 0.992;
be	= 0.1;
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

%% solve for the steady state with the original parameters (just to see that it actually works)
% 
% sol_wcPa_ss([.5,.4])
% %pos_solwcPa = @(logwcPa) sol_wcPa_ss(exp(logwcPa));
% pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);
% 
% 
% [logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)]);
% % loop on tau
% tauH = .1;tauL=0.;
% for itertau = 1:100
% 	tau = 0.5*tauH+0.5*tauL;
% 	[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,logssp,optimset('Display','off'));
% 	wcPa_ss = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
% 	% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
% 	[excess_ss,sseconomy] = sol_wcPa_ss([wcPa_ss]);
% 	budget_def = be*sseconomy(2) - wcPa_ss(1)*tau*(1-sseconomy(2)-sseconomy(1));
% 	if(abs(budget_def)<1e-6 || (tauH-tauL)<1e-6)
% 		break;
% 	elseif (budget_def < 0)
% 		tauH=tau;
% 	elseif(budget_def > 0)
% 		tauL=tau;
% 	end
% end
% 	
%% calibrate it

% initial guesses:
Aa_undevd = .5;
Ym_undevd = 1.;
be_undevd = be/4;
Aa_devd   = 5;
Ym_devd   = 1.;


for cal_iter=1:2%10
	
	%% for developed calibrate abar and Amf so that matches Na_target in ag,
	% Pa_target as the price of agriculture and u_target unemployment by 
	% changing Aa, Ym and Amf

	% to change the calibration target values change Na_target, u_target
	Na_devd_target = 0.1;
	u_devd_target  = 0.07;
	Pa_devd_target = 0.8;
	be = 0.4;
	%cal_devd_fn = @(AaAmfYm) cal_devd_AaYm(AaAmfYm,Na_devd_target,u_devd_target,Pa_devd_target);
	%[x,fval_cal1,exitflag_cal1,out1] = fminsearch(cal_devd_fn,log([Aa_devd, Ym_devd, Amf]));

	cal_devd_fn = @(AaAmfYm) calls_devd_AaYm(AaAmfYm,Na_devd_target,u_devd_target,Pa_devd_target);
	[x,fval_cal1,resid1,exitflag_cal1,out1] = lsqnonlin(cal_devd_fn,[Aa_devd, Ym_devd, Amf],[0,0,0],[10,10,10]);

	
	pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);
	[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)]);
	wcPa_devd = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
	[excess_devd,devd_economy] = sol_wcPa_ss(wcPa_devd);
	devd_logssp = logssp;

	be_devd = be;
	Aa_devd = Aa;
	Ym_devd  = Ym;

	%%

	% for developing calibrate it to Na_target in agriculture, Pa_target as 
	% for the price of agriculture and u_target unemployment by manipulating Aa 
	% Ym and abar
	
	
	% this is the initial guess, for below, when we calibrate to the
	% transition
	Ym	= 1.0;
	Ym_undevd = Ym;

	abar_old = abar;
	
	% to change the calibration target values change Na_target, u_target,
	% Pa_target
	Na_undevd_target = 0.67;
	u_undevd_target  = 0.07;
	Pa_undevd_target = 1.5;
	%cal_undevd_fn = @(abarAaYm) cal_undevd_AaYm(abarAaYm,Na_undevd_target,u_undevd_target,Pa_undevd_target);
	%[x,fval_cal2,exitflag_cal2,out2] = fminsearch(cal_undevd_fn,log([Aa_undevd,be_undevd,abar]));
	cal_undevd_fn = @(abarAaYm) calls_undevd_AaYm(abarAaYm,Na_undevd_target,u_undevd_target,Pa_undevd_target);
	[x,fval_cal2,resid2,exitflag_cal2,out2] = lsqnonlin(cal_undevd_fn,([Aa_undevd,be_undevd,abar]),[0,0,0],[10,Ym,Ym]);
	
	
	pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);
	[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)]);
	wcPa_undevd = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
	[excess_undevd_ss,undevd_ss_economy] = sol_wcPa_ss(wcPa_undevd);

	undevd_ss_logssp = logssp;
	Aa_undevd_ss = Aa;
	be_undevd_ss = be;
	
	if abs(abar -abar_old)<1e-3
		break;
	end
end

%As a first guess, set undevd economy to be the undevd_ss_economy, even
%though it will not be a steady state when we compute the transition

undevd_economy = undevd_ss_economy;
undevd_logssp = undevd_ss_logssp;
p0_trans = [linspace((atan(undevd_logssp(1))+pi/2)*Ym/pi,(atan(devd_logssp(1))+pi/2)*Ym/pi,TT);...
	   linspace(exp(undevd_logssp(2)),exp(devd_logssp (2)),TT)]';


%% transition backwards 
%  First, transition using the approximation that unemployment is the steady state.  
%  This doesn't matter explicitly for any of the choices, just market clearing.

Aa = Aa_devd;

%Here is where you can change the rate of change of the two paths, making
%rate change power over 1 makes for a slow transition and less than 1 is a
%fast transition

rt_chng_pwr	= 1;
Aa_path_chng	= linspace(0,1,TT).^rt_chng_pwr;
Aa_path		= (1-Aa_path_chng)*Aa_undevd_ss+Aa_path_chng*Aa_devd;

rt_chng_pwr	= 1;
Ym_path_chng	= linspace(0,1,TT).^rt_chng_pwr;
Ym_path		= (1-Ym_path_chng)*Ym_undevd+Ym_path_chng*Ym_devd;

be_path = (be_devd-be_undevd)*(Ym_path- Ym_undevd)/(Ym_devd - Ym_undevd) + be_undevd;

trans_path  = zeros(TT,size(devd_economy,2));
excess_path = zeros(TT,2);

% initially hold it to p0
price_path  = p0_trans;
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
for trans_iter =1:1%10
%%
logwA = [devd_logssp(1) log(Aa_devd)];
wcAa_t= [wcPa_devd(1) Aa_devd];
for t = TT-1:-1:1
	Aa = Aa_path(t);
	logwA(2) = log(0.5*Aa + 0.5*exp(logwA(2)));
	logwA(1) = tan( (0.5*price_path(t,1)+ 0.5*wcAa_t(1) )*pi/Ym-pi/2);
	Ym = Ym_path(t);
	be = be_path(t);
	Pa = price_path(t,2);
	pos_solwcAa = @(wcAa) sol_wcAa([(atan(wcAa(1))+pi/2)*Ym/pi exp(wcAa(2))],Pa,trans_path(t+1,:),trans_path(t,2) );
	tauH = 0.05; tauL=0.001;
	for itertau = 1:100
		tau = 0.5*tauH+0.5*tauL;
		p0 = logwA;
		[logwA, fval,exitflag,output,J] = fsolve(pos_solwcAa,p0,optimset('Display','off'));
		wcAa_t = [(atan(logwA(1))+pi/2)*Ym/pi exp(logwA(2))];
		% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
		[excess_trans,trans_economy] = sol_wcAa(wcAa_t,Pa,trans_path(t+1,:),trans_path(t,2));
		budget_def = be*trans_economy(2) - wcAa_t(1)*tau*(1-trans_economy(2)-trans_economy(1));
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
	solpath_back(t,:) = logwA;
	price_path_back(t,:) = [wcAa_t(1) Pa];
	Aa_implied_back(t) = wcAa_t(2);
end
price_path(:,1) = price_path_back(:,1); % replace wages, but hold fixed the Pa;
Aa_implied_fwd(1) = wcAa_t(2);

%% now go forward to get quantities right
for t = 2:TT-1
	Aa = Aa_path(t);
	Ym = Ym_path(t);
	be = be_path(t);
	Pa = price_path(t,2);
	if t>2
		pos_solwcAa = @(wcAa) sol_wcAa_fwd([(atan(wcAa(1))+pi/2)*Ym/pi exp(wcAa(2))],Pa,trans_path(t+1,:),trans_path(t-1,1:2));
	else
		pos_solwcAa = @(wcAa) sol_wcAa_fwd([(atan(wcAa(1))+pi/2)*Ym/pi exp(wcAa(2))],Pa,trans_path(t+1,:),[Na_undevd_target,u_undevd_target*(1-Na_undevd_target)]);
	end
	tauH = 0.05; tauL=0.001;
	for itertau = 1:100
		tau = 0.5*tauH+0.5*tauL;

		p0 = solpath_back(t,:); 
		[logwA, fval,exitflag,output,J] = fsolve(pos_solwcAa,p0,optimset('Display','off'));
		wcAa_t = [(atan(logwA(1))+pi/2)*Ym/pi exp(logwA(2))];
		% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
		if(exitflag<=0)
			logwA = p0;
			wcAa_t = [price_path(t,1)  Aa_path(t)];
		end
		
		if t>2
			[excess_trans,trans_economy] = sol_wcAa_fwd(wcAa_t,Pa,trans_path(t+1,:),trans_path(t-1,1:2));
		else
			[excess_trans,trans_economy] = sol_wcAa_fwd(wcAa_t,Pa,trans_path(t+1,:),[Na_undevd_target,(1-Na_undevd_target)*u_undevd_target]);
		end
		budget_def = be*trans_economy(2) - wcAa_t(1)*tau*(1-trans_economy(2)-trans_economy(1));
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
	price_path_fwd(t,:) = [wcAa_t(1) Pa]; % or if need to go slower: wcPa_t*.25 + .75*price_path_back(t,:);
	Aa_implied_fwd(t) = wcAa_t(2);
end
trans_path(TT,:) = devd_economy;
price_dif = abs(price_path_fwd - price_path_back);

% check the calibration criteria at our calibration period.  
	theeconomy	= trans_path(1,:);
	calresid(1)	= Na_undevd_target - theeconomy(1);
	calresid(2)	= u_undevd_target - theeconomy(2)/(1-theeconomy(1));
	%this makes prices the target.  
	%calresid(3)	= Pa_undevd_target - price_path(1,2);
	%Instead I impose prices, and make sure the Aa fits:
	calresid(3)	= Aa_path(1) - Aa_implied_fwd(1);

	resid = calresid.^2;
	
% adjust the path of ag TFP growth if price path is off:
abs(price_path_back(:,2) - p0_trans(:,2));

end




%% Compute paths for revenue per worker in Ag & Urban at P_t and P_0

rev_pt = [price_path(:,2).*Aa.*trans_path(:,1).^mu Ym*ones(TT,1)];  % ag,man
percaprev_pt = [price_path(:,2).*Aa.*trans_path(:,1).^(mu-1) Ym./(1-trans_path(:,1))];  % ag,man
rev_p0 = [price_path(1,2).*Aa*trans_path(:,1).^mu Ym*ones(TT,1)];  % ag,man
percaprev_p0 = [price_path(1,2).*Aa*trans_path(:,1).^(mu-1) Ym./(1-trans_path(:,1))];  % ag,man
rev_pTT = [price_path(TT,2).*Aa*trans_path(:,1).^mu Ym*ones(TT,1)];  % ag,man
percaprev_pTT = [price_path(TT,2).*Aa*trans_path(:,1).^(mu-1) Ym./(1-trans_path(:,1))];  % ag,man


%%
save_plots =1;
cd trans_AaYm_results/linAa_veryslowYm

save trans_space
upath = trans_path(:,2)./(1-trans_path(:,1));
mpath = (-trans_path(2:TT,1) + trans_path(1:TT-1,1))./trans_path(1:TT-1,1);
%
figure(1);
[ax,h1,h2]=plotyy([1:TT],trans_path(:,1),[1:TT],Aa_path);title('Fraction in agriculture');
set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','N_a','A_a');
if (save_plots == 1) saveas(gca,'Natrans','eps2c'); end

figure(2);
[ax,h1,h2]=plotyy([2:TT],mpath,[1:TT],Aa_path);title('Rural -> urban rate');
set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','m/N_a','A_a');
if (save_plots == 1) saveas(gca,'mtrans','eps2c'); end

figure(3);
[ax,h1,h2]=plotyy([1:TT],upath,[1:TT],Ym_path);title('Unemployment rate');
set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','u','Y_m');
if (save_plots == 1) saveas(gca,'utrans','eps2c'); end

figure(4);
[ax,h1,h2]=plotyy([1:TT],Amf*trans_path(:,3).^(1-eta),[1:TT],Ym_path);title('Job finding rate');
set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','p(Q)','Y_m');
if (save_plots == 1) saveas(gca,'pQtrans','eps2c'); end
	
figure(5);
[ax,h1,h2]=plotyy([1:TT],percaprev_pt(:,1)./percaprev_pt(:,2) ,[1:TT], percaprev_pt(:,1));title('Relative Revenue Per Capita, P_t');
set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','NorthEast','P_t y_a/y_m','P_t y_a');
y20 = min(percaprev_pt(:,1));y2Y=max(percaprev_pt(:,1));
y10 = min(percaprev_pt(:,1)./percaprev_pt(:,2));y1Y=max(percaprev_pt(:,1)./percaprev_pt(:,2));
set(ax(1),'YTick',round(([0:6]*(y1Y-y10)/6 +y10)*100)/100 );set(ax(2),'YTick',round(([0:6]*(y2Y-y20)/6 +y20)*100)/100);
if (save_plots == 1) saveas(gca,'rev_Pt_trans','eps2c'); end

figure(6);
h=plot([1:TT],percaprev_p0(:,1)./percaprev_p0(:,2) ,[1:TT], percaprev_p0(:,1));title('Relative Revenue Per Capita, P_0');
set(h,'LineWidth',2);legend('Location','North','P_0 y_a/y_m','P_0 y_a');
if (save_plots == 1) saveas(gca,'rev_P0_trans','eps2c'); end


figure(7);
h=plot([1:TT],percaprev_pTT(:,1)./percaprev_pTT(:,2) ,[1:TT], percaprev_pTT(:,1));title('Relative Revenue Per Capita, P_T');
set(h,'LineWidth',2);legend('Location','North','P_T y_a/y_m','P_T y_a');
if (save_plots == 1) saveas(gca,'rev_PTT_trans','eps2c'); end

figure(8);
[ax,h1,h2]=plotyy([1:TT],price_path(:,2),[1:TT],Aa_path);title('Relative price of agricultural good');
set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','P_t','A_a');
if (save_plots == 1) saveas(gca,'agPrice','eps2c'); end


cd ~/Documents/CurrResearch/Devt/Computation
