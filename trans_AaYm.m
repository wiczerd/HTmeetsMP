% solves transition paths from low to high steady states

% This program computes the transition for Harris-Todaro Meets
% Mortensen-Pissarides

cd ~/Documents/CurrResearch/Devt/Computation

global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau

TT_data = 400;
TT_leadin = 100;
TT = TT_data+TT_leadin;
save_plots =0;

cbar	= 0;%-0.6; % note this is the inverse because I changed the util function
abar	= 0.09; % this gets changed below in the calibration
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

% initial guesses:
Aa_undevd = .5;
Ym_undevd = 1.;
be_undevd = be/4;
Aa_devd   = 5;
Ym_devd   = 1.;


for cal_iter=1:4
	
	%% for developed calibrate abar and Amf so that matches Na_target in ag,
	% Pa_target as the price of agriculture and u_target unemployment by 
	% changing Aa, Ym and Amf

	% to change the calibration target values change Na_target, u_target
	Na_devd_target = 0.05;
	u_devd_target  = 0.075;
	Pa_devd_target = 1.0;
	be = 0.4;
	cal_devd_fn = @(AaAmfYm) cal_devd_AaYm(AaAmfYm,Na_devd_target,u_devd_target,Pa_devd_target);
	[x,fval_cal1,exitflag_cal1,out1] = fminsearch(cal_devd_fn,log([Aa_devd, Ym_devd, Amf]));
	%cal_devd_fn = @(AaAmf) cal_devd_AaYm2(AaAmf,Na_target,u_target);
	%[x,fval_cal1,exitflag_cal1,out1] = fminsearch(cal_devd_fn,log([Aa_devd, abar]));

	pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);
	[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)]);
	wcPa_ss = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
	[excess_devd,devd_economy] = sol_wcPa_ss(wcPa_ss);
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
	Na_undevd_target = 0.75;
	u_undevd_target  = 0.07;
	Pa_undevd_target = 1.5;
	cal_undevd_fn = @(abarAaYm) cal_undevd_AaYm(abarAaYm,Na_undevd_target,u_undevd_target,Pa_undevd_target);
	[x,fval_cal2,exitflag_cal2,out2] = fminsearch(cal_undevd_fn,log([Aa_undevd,be_undevd,abar]));
	
	
	pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);
	[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)]);
	wcPa_ss = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
	[excess_undevd_ss,undevd_ss_economy] = sol_wcPa_ss(wcPa_ss);

	undevd_ss_logssp = logssp;
	Aa_undevd_ss = Aa;
	be_undevd_ss = be;
	
	
	if abs(abar -abar_old)<1e-3
		break;
	end
end

%%

% compute a low steady state 100 periods before the calibration period
Na_ss0_target = (Na_devd_target-Na_undevd_target)/TT_data*(-TT_leadin)+Na_undevd_target;
cal_ss0_fn = @(Aahere) cal_ss0_AaYm(Aahere,Na_ss0_target);
[x,fval_cal3,exitflag_cal3,out3] = fminsearch(cal_ss0_fn,log(Aa_undevd));

pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);
[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)]);
wcPa_ss = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
[excess_undevd_ss0,undevd_ss0_economy] = sol_wcPa_ss(wcPa_ss);

ss0_logssp = logssp;
Aa_ss0 = Aa;
Ym_ss0 = Ym;



%% transition backwards 
%  First, transition using the approximation that unemployment is the steady state.  
%  This doesn't matter explicitly for any of the choices, just market clearing.

Aa = Aa_devd;

%Here is where you can change the rate of change of the two paths, making
%rate change power over 1 makes for a slow transition and less than 1 is a
%fast transition

rt_chng_pwr	= 1;
Aa_path_chng	= linspace(0,1,TT).^rt_chng_pwr;
Aa_path		= (1-Aa_path_chng)*Aa_ss0+Aa_path_chng*Aa_devd;

rt_chng_pwr	= 5;
Ym_path_chng	= linspace(0,1,TT).^rt_chng_pwr;
Ym_path		= (1-Ym_path_chng)*Ym_ss0+Ym_path_chng*Ym_devd;

trans_path  = zeros(TT,size(devd_economy,2));
excess_path = zeros(TT,2);
price_path  = zeros(TT,2);


p0_trans = [linspace(exp(devd_logssp(2)),exp(ss0_logssp(2)),TT);...
	   linspace((atan(devd_logssp(1))+pi/2)*Ym/pi,(atan(ss0_logssp(1))+pi/2)*Ym/pi,TT)]';

trans_economy	= devd_economy;
trans_path(TT,:)= trans_economy;

trans_economy	= undevd_ss0_economy;
trans_path(1,:)= trans_economy;


% use steady state unemployment first
trans_path(:,2) = -1.0;

price_path(TT,:)= p0_trans(TT,:);

logssp = devd_logssp;

for trans_iter =1:10

for t = TT-1:-1:1
	Aa = Aa_path(t);
	Ym = Ym_path(t);
	be = (be_devd-be_undevd)*Ym/(Ym_devd - Ym_undevd) + be_undevd;
	pos_solwcPa = @(wcPa) sol_wcPa([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))],trans_path(t+1,:),trans_path(t,2));
	tauH = 0.05; tauL=0.001;
	for itertau = 1:100
		tau = 0.5*tauH+0.5*tauL;
		%p0 = [tan(p0_trans(t,1)*pi/Ym-pi/2) log(p0_trans(t,2))];
		%if(t>TT/2)
		%	p0 = devd_logssp;
		%else
		%	p0 = undevd_logssp;
		%end
		p0 = logssp;
		[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,p0,optimset('Display','off'));
		wcPa_t = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
		% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
		[excess_trans,trans_economy] = sol_wcPa(wcPa_t,trans_path(t+1,:),trans_path(t,2));
		budget_def = be*trans_economy(2) - wcPa_t(1)*tau*(1-trans_economy(2)-trans_economy(1));
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
	price_path(t,:) = wcPa_t;
end

price_path_back = price_path;

% now go forward to get quantities right
for t = 2:TT-1
	Aa = Aa_path(t);
	Ym = Ym_path(t);
	pos_solwcPa = @(wcPa) sol_wcPa_fwd([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))],trans_path(t+1,:),trans_path(t-1,1:2));
	tauH = 0.05; tauL=0.001;
	for itertau = 1:100
		tau = 0.5*tauH+0.5*tauL;
		%p0 = [tan(p0_trans(t,1)*pi/Ym-pi/2) log(p0_trans(t,2))];
% 		if(t>TT/2)
% 			p0 = devd_logssp;
% 		else
% 			p0 = undevd_logssp;
% 		end
		p0 = [tan(price_path_back(t,1)*pi/Ym-pi/2) log(price_path_back(t,2))]; 
		[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,p0,optimset('Display','off'));
		wcPa_t = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
		% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
		[excess_trans,trans_economy] = sol_wcPa_fwd(wcPa_t,trans_path(t+1,:),trans_path(t-1,1:2));
		budget_def = be*trans_economy(2) - wcPa_t(1)*tau*(1-trans_economy(2)-trans_economy(1));
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
	price_path(t,:) = wcPa_t; % or if need to go slower: wcPa_t*.25 + .75*price_path_back(t,:);
end
trans_path(TT,:) = devd_economy;
price_dif = abs(price_path - price_path_back);
if(max(max(price_dif))<1e-6)
	break;
else % this is a heuristic in case it's not converging
	disp('Max deviation:')
	max(max(price_dif))
	for fixiter=1:10
		papprox_coef = quantreg([1:TT]',price_path(:,2),.5,5,100);
		%plot([1:TT],price_path(:,2),[1:TT],polyval(papprox_coef,[1:TT]'));

		pbad = zeros(TT,1)==1;
		% replace values that deviate from the 5th-order fit and have discontinuous
		% derivatives
		papprox = polyval(papprox_coef,time);
		for t=2:TT-1
			if(abs(papprox(t)-price_path(t,2))/(1+papprox(t))>.05 && abs(price_path(t+1,2)-2*price_path(t,2)+price_path(t-1,2))>.05 ) 
				pbad(t)=1;
			end
			if(trans_path(t,1)<.01)
				pbad(t)=1;
			end
		end
		pgood = pbad==0;
		if(sum(pbad)>0)
			price_path(pbad,2) = interp1(time(pgood),price_path(pgood,2),time(pbad));
		else
			break;
		end
	end
end

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
