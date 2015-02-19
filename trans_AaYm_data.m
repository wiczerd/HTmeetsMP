% solves transition paths from low equilibrium to high steady state

% This program computes the transition for Harris-Todaro Meets
% Mortensen-Pissarides

cd ~/Documents/CurrResearch/Devt/Computation

global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau

TT = 150*4;
save_plots =1; %save the plots?
data_plots =0; %plot the data?
param_update = .75;

optsoff = optimset('Display','off');

cbar	= 0;%-0.6; % note this is the inverse because I changed the util function
abar	= 0.4; % this gets changed below in the calibration
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
be	= 0.4;
tau	= 0.0;

%import the price sequence from Alvarez-Cuadrado & Poschke:

ts_data = dir('*_agprod.csv');
Ncountry = size(ts_data,1);
c_name = cell(Ncountry,1);
c_nqtr = zeros(Ncountry,1);
c_Pa = zeros(151,Ncountry);
c_time = zeros(151,Ncountry);
c_Pa_qtr = zeros(150*4,Ncountry);
c_time_qtr= zeros(150*4,Ncountry);

for ci = 1:Ncountry
	c_name(ci) = cellstr(strtok(ts_data(ci).name,'_'));
	fid = fopen(ts_data(ci).name);
	for li = 1:151
		c_agprod_l = fgetl(fid);
		mis = strsplit(c_agprod_l,' ');
		mis0= sscanf(char(mis(end)),'%d');
		mis1= sscanf(char(mis(end-1)),'%f');
		mis2= sscanf(char(mis(end-2)),'%f');
		if(~isempty(mis1) && ~isempty(mis2))
			
			c_Pa(li,ci) = 1/mis2;
			c_time(li,ci) = mis0;
		end
	end
	indic_data = c_time(:,ci)>0;
	nyr = sum(indic_data);
	yr0 = min(c_time(indic_data,ci));
	Pa_qtr = interp1(c_time(indic_data,ci)-yr0,...
						c_Pa(indic_data,ci), 0:0.25:nyr-1);
	Pa_qtr = Pa_qtr(1:end-1);
	c_nqtr(ci) = size(Pa_qtr,2);
	c_Pa_qtr(1:c_nqtr(ci),ci) = Pa_qtr';
	time_qtr = interp1(c_time(indic_data,ci)-yr0,...
						c_time(indic_data,ci), 0:0.25:nyr-1);
	c_time_qtr(1:c_nqtr(ci),ci) = time_qtr(1:end-1)';
	
	if(data_plots ==1)
		h=plot(c_time(indic_data,ci),c_Pa(indic_data,ci) );
		title(['Relative Price of Agricultural Goods, ' c_name{ci}],'FontSize',14);
		set(h,'LineWidth',2);
		set(gcf,'color','white');
		grid on;
		if save_plots==1
			saveas(gca,['relprice_' c_name{ci} '.eps'],'eps2c');
			saveas(gca,['relprice_' c_name{ci} '.png']); 
		end
	end
end

% AC-P calibration targets data: 
cal_period = ones(Ncountry,1);
cal_Na0 = ones(Ncountry,1);
cal_NaTT = ones(Ncountry,1);
cal_uTT  = ones(Ncountry,1);
% Canada
cal_period(1) = 1;
cal_Na0(1)	= 0.458;
cal_NaTT(1)	= 0.1985*1/6 + 0.157*5/6; % linear interpolation on AC-P data.
cal_uTT(1)	= 0.028; %(2.3 + 2.8 + 3.6 + 2.4 + 2.9)/5
% Germany
cal_period(2) = find(c_time_qtr(:,2)==1849);
cal_Na0(2)	= 0.5601836;
cal_NaTT(2)	= 0.3455502;
cal_uTT(2)	= 0.0504; %(4.9+4.3*2+5.7+6.0)/5 --- 5 year moving average
% UK
cal_period(3) = 1;
cal_Na0(3)	= 0.37;
cal_NaTT(3)	= 0.1;
cal_uTT(3)	= 0.035; %(3.1+3.2+4.2)/3 -- 3 year moving average centered on 1912
% US
cal_period(4) = 1;
cal_Na0(4)	= 0.73;
cal_NaTT(4)	= 0.121;
cal_uTT(4)	= 0.0426; %(3.8+5.9+5.3+3.3+3.0)/5



%Gollin data
apg_gol = zeros(7,2,4);
% Canada
apg_gol(:,1,1) = [1925 1930 1940 1950 1960 1970 1980];
apg_gol(:,2,1) = [0.401514604 0.2843524543 0.400050586 0.579004329 0.462501694 0.7100280823 0.5554019619];
% Germany
apg_gol(:,1,2) = [1880 1895 1905 1925 1930 1935 1950];
apg_gol(:,2,2) = [0.640965982 0.676163922 0.6023711007 0.4340359094 0.6120689655 0.5435361042 0.4824311491];
%UK
apg_gol(:,1,3) = [1840 1850 1860 1870 1880 1890 1900];
apg_gol(:,2,3) = [0.9838880449 0.9413651571 0.9456268448 0.9776863002 0.8134439608 0.8228239608 0.7137124952];
%US
apg_gol(:,1,4) = [1869 1879 1899 1907 1918 1926 1950];
apg_gol(:,2,4) = [0.3054335851 0.2451210583 0.3804705773 0.5510643855 0.7142943069 0.4829975482 0.5763560775];

for ci=1:Ncountry
	indic_data = c_time(:,ci)>0;
	if(data_plots ==1)
		h=plot(c_time(indic_data,ci),c_Pa(indic_data,ci),apg_gol(:,1,ci),apg_gol(:,2,ci));
		title(['APG,' c_name{ci}],'FontSize',14);legend('Location','SouthWest','P_t','APG');
		set(h,'LineWidth',2);
		set(gcf,'color','white');
		grid on;
		if save_plots==1
			saveas(gca,['relprice_golprod_' c_name{ci} '.eps'],'eps2c');
			saveas(gca,['relprice_golprod_' c_name{ci} '.png']); 
		end
	end
end
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

 	
%% calibrate it for country ci

%for ci = 1:Ncountry
ci = 2;
	TT = c_nqtr(ci);

	% initial guesses:
	Aa_undevd_ss = .5;
	Ym_undevd = 1.;
	be_undevd_ss = be;
	Aa_devd   = 5;
	Ym_devd   = 1.;
	Amf_devd  = Amf;

	for cal_iter=1:20

		%% for developed calibrate abar and Amf so that matches Na_target in ag,
		% Pa_target as the price of agriculture and u_target unemployment by 
		% changing Aa, Ym and Amf

		% to change the calibration target values change Na_target, u_target
		Na_devd_target = cal_NaTT(ci);
		u_devd_target = cal_uTT(ci);
		Pa_devd_target = mean(c_Pa_qtr(TT-20:TT,ci));

		Amf_old = Amf;

		cal_devd_fn = @(AaAmfYm) calls_devd_AaYm(AaAmfYm,Na_devd_target,u_devd_target,Pa_devd_target);
		[x,fval_cal1,resid1,exitflag_cal1,out1] = lsqnonlin(cal_devd_fn,[Aa_devd, Ym_devd, Amf_devd],[0,0,0],[20,20,10],optsoff);
		if exitflag_cal1<1
			disp(exitflag_cal1)
			disp(fval_cal1)
		end
			
		Amf = param_update*Amf +(1-param_update)*Amf_old; 

		pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);
		[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)],optsoff);
		wcPa_devd = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
		[excess_devd,devd_economy] = sol_wcPa_ss(wcPa_devd);
		devd_logssp = logssp;

		if exitflag<1
			disp(exitflag);
			disp(fval);
		end
		
		be_devd = be;
		Aa_devd = Aa;
		Ym_devd  = Ym;
		Amf_devd = Amf;
		
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
		Na_undevd_target = cal_Na0(ci);
		Pa_undevd_target = mean(c_Pa_qtr(1:20,ci));
		
		if cal_period(ci) == 1
			Na_undevd_target_ss = Na_undevd_target;
		else
				% half way between assuming linear extrapolation and
				% assuming no change
			Na_undevd_target_ss = (cal_Na0(ci)-(cal_NaTT(ci)-cal_Na0(ci))/(c_nqtr(ci)- cal_period(ci))*cal_period(ci))*0.5 + .5*Na_undevd_target;
		end
		
	%	cal_undevd_fn = @(abarAaYm) calls_undevd_AaYm(abarAaYm,Na_undevd_target,u_undevd_target,Pa_undevd_target);
	%	[x,fval_cal2,resid2,exitflag_cal2,out2] = lsqnonlin(cal_undevd_fn,[Aa_undevd,be_undevd_ss,abar],[0,0,0],[10,Ym,Ym]);
		cal_undevd_fn = @(abarAaYm) calls_undevd_AaYm_fixbe(abarAaYm,Na_undevd_target_ss,Pa_undevd_target);
		[x,fval_cal2,resid2,exitflag_cal2,out2,lam2, jac2] = lsqnonlin(cal_undevd_fn,[Aa_undevd_ss,abar],[0,0],[10,Ym],optsoff);
		if exitflag_cal2<1
			disp(exitflag_cal2)
			disp(fval_cal2)
		end
		abar = abar*param_update + (1-param_update)*abar_old;
		
		pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);
		[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)],optsoff);
		wcPa_undevd = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
		[excess_undevd_ss,undevd_ss_economy] = sol_wcPa_ss(wcPa_undevd);

		if exitflag<1
			disp(exitflag);
			disp(fval);
		end
		
		undevd_ss_logssp = logssp;
		Aa_undevd_ss = Aa;
		Aa_undevd = Aa_undevd_ss;
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
	p0_trans(1:TT,2) = c_Pa_qtr(1:TT,ci);

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

	%rt_chng_pwr	= 1;
	%Ym_path_chng	= linspace(0,1,TT).^rt_chng_pwr;
	%Ym_path	= (1-Ym_path_chng)*Ym_undevd+Ym_path_chng*Ym_devd;
	rt_chng_pwr	= (Ym_devd/Ym_undevd)^(1/(TT-1));
	Ym_path		= Ym_undevd*rt_chng_pwr.^(linspace(0,TT-1,TT));

	be_path = (be_devd-be_undevd_ss)*(Ym_path- Ym_undevd)/(Ym_devd - Ym_undevd) + be_undevd_ss;

	trans_path  = zeros(TT,size(devd_economy,2));
	excess_path = zeros(TT,2);
	uss_path = zeros(TT,1);
	utback_path= zeros(TT,1);
	tau_path = zeros(TT,1);
	bad_periods = zeros(TT,1);

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
	solpath_fwd= zeros(size(price_path));
	Aa_implied_back = ones(TT,1)*Aa_devd;
	Aa_implied_fwd  = ones(TT,1)*Aa_devd;


	price_path_back(TT,:) = wcPa_devd;
	price_path_fwd(TT,:) = wcPa_devd;
	for trans_iter =1:100
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
				if exitflag<1 %try with a new starting point
					p0 = solpath_back(t+1,:);
					[logwA, fval,exitflag,output,J] = fsolve(pos_solwcAa,p0,optimset('Display','off'));
				end
				
				wcAa_t = [(atan(logwA(1))+pi/2)*Ym/pi exp(logwA(2))];
				% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
				[excess_trans,trans_economy] = sol_wcAa(wcAa_t,Pa,trans_path(t+1,:),trans_path(t,2));
				%absolute level of UI replacement
				%budget_def = be*theeconomy(2) - wcPa_ss(1)*tau*(1-theeconomy(2)-theeconomy(1));
				%propotional UI replacement
				budget_def = be*trans_economy(2) - tau*(1-trans_economy(2)-trans_economy(1));
				if(abs(budget_def)<1e-6 || (tauH-tauL)<1e-6)
					break;
				elseif (budget_def < 0)
					tauH=tau;
				elseif(budget_def > 0)
					tauL=tau;
				end
			end
			if trans_economy(1) >= 1. || excess_trans*excess_trans' >1e-4;
				trans_economy = trans_path_back(t+1,:);
				bad_periods(t) = 1;
				solpath_back(t,:) = solpath_back(t+1,:);
			else
				% this is the guess for the path forward
				solpath_back(t,:) = logwA;
			end
			tau_path(t) = tau;
			trans_path_back(t,:) = trans_economy;
			trans_path(t,:) = trans_economy;
			excess_path(t,:)= excess_trans;
			
			price_path_back(t,:) = [wcAa_t(1) Pa];
			Aa_implied_back(t) = wcAa_t(2);
			pQ	= Amf*trans_economy(3)^(1-eta);
			uss_path(t) = lambda*(1-trans_economy(1))*(1-pQ)/(pQ + lambda*(1-pQ) );
		end
		price_path(:,1) = price_path_back(:,1); % replace wages, but hold fixed the Pa;
		Aa_implied_fwd(1) = wcAa_t(2);
		price_path_fwd(1,:) = [wcAa_t(1) Pa];
		utback_path(1) = trans_path_back(1,2)*(1-trans_path_back(1,1));
		for t = 2:TT-1
			pQ	= Amf*trans_path_back(t,3)^(1-eta);
			Na_t = trans_path_back(t,1);
			Na_tm1 = trans_path_back(t-1,1);
			u_tm1 = utback_path(t-1);
			utback_path(t) = (1-pQ)*(u_tm1 + lambda*(1-u_tm1 -Na_t) + (Na_tm1 - Na_t));
		end
	%	trans_path(:,2) = utback_path;
		%% now go forward to get quantities right
		for t = 2:TT-1
			Aa = Aa_path(t);
			Ym = Ym_path(t);
			be = be_path(t);
			Pa = price_path(t,2);

			pos_solwcAa = @(wcAa) sol_wcAa_fwd([(atan(wcAa(1))+pi/2)*Ym/pi exp(wcAa(2))],Pa,trans_path(t+1,:),trans_path(t-1,1:2));

			tau = tau_path(t);
		
			p0 = solpath_back(t,:); 
			[logwA, fval,exitflag,output,J] = fsolve(pos_solwcAa,p0,optimset('Display','off'));
			wcAa_t = [(atan(logwA(1))+pi/2)*Ym/pi exp(logwA(2))];
			if(exitflag<=0) % try again with a different starting point
				p0 = solpath_fwd(t-1,:); 
				[logwA, fval,exitflag,output,J] = fsolve(pos_solwcAa,p0,optimset('Display','off'));
				wcAa_t = [(atan(logwA(1))+pi/2)*Ym/pi exp(logwA(2))];
			end
			if(exitflag<=0) %it still didn't work
				logwA = p0;
				wcAa_t = [price_path_fwd(t-1,1)  Aa_implied_fwd(t-1)];
			end
			solpath_fwd(t,:) = logwA;
			% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
			[excess_trans,trans_economy] = sol_wcAa_fwd(wcAa_t,Pa,trans_path(t+1,:),trans_path(t-1,1:2));
			budget_def = be*trans_economy(2) - tau*(1-trans_economy(2)-trans_economy(1));
			if trans_economy(1) >= 1. || exitflag<0;
				trans_economy = trans_path(t-1,:);
				bad_periods(t) = 1;
			end
			
			trans_path(t,:) = trans_economy;
			excess_path(t,:)= excess_trans;
			price_path_fwd(t,:) = [wcAa_t(1) Pa]; % or if need to go slower: wcPa_t*.25 + .75*price_path_back(t,:);
			Aa_implied_fwd(t) = wcAa_t(2);
		end
		trans_path(TT,:) = devd_economy;
		price_dif = abs(price_path_fwd - price_path_back);
		price_path = price_path_fwd;
	% check the calibration criteria at our calibration period.  

		abar_old = abar;
		be_old = be_path;
	%%	
		theeconomy	= trans_path(cal_period(ci),:);
		calresid(1)	= Na_undevd_target - theeconomy(1);
		% this makes the unemployment target
		%calresid(2)	= u_undevd_target - theeconomy(2)/(1-theeconomy(1));
		%this makes prices the target.  
		calresid(2)	= Pa_undevd_target - price_path(1,2);

		%Convergence in A? :
		calresid(3)	= sum(abs(Aa_implied_back - Aa_implied_fwd));

		resid = calresid.^2;

		Aa_path = param_update*Aa_implied_fwd'+(1-param_update)*Aa_path;

		%re-calibrate the low-development state with just abar
		Aa = Aa_path(1);
		Ym = Ym_path(1);
		Pa = price_path(1,2);
		w0 = price_path(1,1);
		be = be_path(1);
		options = optimoptions('fsolve','Jacobian','off');
		cal_undevd_fn = @(Aaabar) calls_undevd_AaYm_bkwd(Aaabar,Pa,w0,Na_undevd_target,trans_path(cal_period(ci)+1,:),trans_path(cal_period(ci),2));
		[x,fval_caliter,resid_caliter,exitflag_caliter,out_caliter,J_caliter] = lsqnonlin(cal_undevd_fn,[Pa*1.1 abar*.9],[0 0],[10 Ym]);
		%[x,fval_caliter,resid_caliter,exitflag_caliter,out_caliter,J_caliter] = lsqnonlin(cal_undevd_fn,abar,0,Ym);

		%pos_solwc = @(wc) sol_wc((atan(wc)+pi/2)*Ym/pi,Pa,trans_path(2,:),trans_path(1,2)); % should I use trans_path(1,2)) instead of utarget*(1-Natarget)
		pos_solwc = @(wcPa) sol_wcPa([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))],trans_path(2,:),trans_path(1,2)); % should I use trans_path(1,2)) instead of utarget*(1-Natarget)
		[tanw, fval,exitflag,output,J] = fsolve(pos_solwc,[tan(w0*pi/Ym-pi/2) log(Pa)]);
		[excess_undevd,undevd_economy] = sol_wcPa([(atan(tanw(1))+pi/2)*Ym/pi exp(tanw(2))],trans_path(2,:),trans_path(1,2));
		trans_path(1,:) = undevd_economy;

	 %%

		param_resid = abs(be-be_path(1)) + abs(abar - abar_old);
		be_undevd = be;
	%	be_undevd = be_path(1);
		be_path_implied = (be_devd-be_undevd)*(Ym_path- Ym_undevd)/(Ym_devd - Ym_undevd) + be_undevd;
		be_path = param_update*be_path_implied + (1-param_update)*be_path;
		abar = param_update*abar + (1-param_update)*abar_old;

		if (sum(resid)<1e-6) || (trans_iter >2 && param_resid< 1e-6)
			break;
		end
		bad_periods(:) = 0;

	end




	%% Compute paths for revenue per worker in Ag & Urban at P_t and P_0

	rev_pt = [price_path(:,2).*Aa_path'.*trans_path(:,1).^mu Ym_path'];  % ag,man
	percaprev_pt = [price_path(:,2).*Aa_path'.*trans_path(:,1).^(mu-1) Ym_path'./(1-trans_path(:,1)-trans_path(:,2))];  % ag,man
	rev_p0 = [price_path(1,2).*Aa_path'.*trans_path(:,1).^mu Ym_path'];  % ag,man
	percaprev_p0 = [price_path(1,2).*Aa_path'.*trans_path(:,1).^(mu-1) Ym_path'./(1-trans_path(:,1)-trans_path(:,2))];  % ag,man
	% this is using the TT period prices: fully developed
	rev_pTT = [price_path(TT,2).*Aa_path'.*trans_path(:,1).^mu Ym_path'];  % ag,man
	percaprev_pTT = [price_path(TT,2).*Aa_path'.*trans_path(:,1).^(mu-1) Ym_path'./(1-trans_path(:,1)-trans_path(:,2) )];  % ag,man


	%%
	mkdir('trans_AaYm_results',['calAa_' c_name{ci} 'Pa_linYm']);
	cd(['trans_AaYm_results/calAa_' c_name{ci} 'Pa_linYm'])

	save(['trans_space_' c_name{ci} '.mat']);

	%%
	upath = trans_path(:,2)./(1-trans_path(:,1));
	upath_back = trans_path_back(:,2)./(1-trans_path_back(:,1));

	mpath = (-trans_path(2:TT,1) + trans_path(1:TT-1,1))./trans_path(1:TT-1,1);
	%
	figure(1);
	[ax,h1,h2]=plotyy([1:TT],trans_path(:,1),[1:TT],Aa_path);title('Fraction in agriculture','FontSize',14);
	set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','N_a','A_a');
	set(gcf,'color','white');
	set(ax(1), 'YLim', [0.0 1.0]);
	set(ax(1), 'YTick', [0.0:0.2:1.0]);
	set(ax(2), 'YLim', [1.0 7.]);
	set(ax(2), 'YTick', [1.0:1.5:7.]);
	grid on;
	if (save_plots == 1) 
		saveas(gca,'Natrans','eps2c');
		saveas(gca,'Natrans.png'); 
	end

	figure(2);
	[ax,h1,h2]=plotyy([2:TT],mpath,[1:TT],Aa_path);title('Rural -> urban rate','FontSize',14);
	set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','m/N_a','A_a');
	set(gcf,'color','white');
	%axes(ax(1));axis([-inf inf 0.0 0.001]);
	set(ax(1), 'YLim', [0.0 0.02]);
	set(ax(1), 'YTick', [0.0:0.002:0.02]);
	grid on;
	if (save_plots == 1) 
		saveas(gca,'mtrans','eps2c'); 
		saveas(gca,'mtrans.png'); 
	end

	figure(3);
	[ax,h1,h2]=plotyy([1:TT],upath,[1:TT],Ym_path);title('Unemployment rate','FontSize',14);
	set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','u','Y_m');
	set(gcf,'color','white');
	set(ax(1), 'YLim', [0.05 0.15]);
	set(ax(1), 'YTick', [0.05:0.025:0.15]);
	set(ax(2), 'YLim', [0.0 10.0]);
	set(ax(2), 'YTick', [0.0:2.5:10.0]);
	grid on;
	if (save_plots == 1) 
		saveas(gca,'utrans','eps2c'); 
		saveas(gca,'utrans.png'); 
	end

	figure(4);
	[ax,h1,h2]=plotyy([1:TT],Amf*trans_path(:,3).^(1-eta),[1:TT],Ym_path);title('Job finding rate','FontSize',14);
	set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','p(Q)','Y_m');
	set(gcf,'color','white');
	grid on;
	if (save_plots == 1) 
		saveas(gca,'pQtrans','eps2c'); 
		saveas(gca,'pQtrans.png'); 
	end

	figure(5);
	[ax,h1,h2]=plotyy([2:TT],percaprev_pt(2:end,1)./percaprev_pt(2:end,2) ,[1:TT], percaprev_pt(:,1));title('Relative Revenue Per Capita, P_t','FontSize',14);
	set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','South','P_t y_a/y_m','P_t y_a');
	y20 = min(percaprev_pt(:,1));
	y2Y=max(percaprev_pt(:,1));
	y10 = min(percaprev_pt(:,1)./percaprev_pt(:,2));
	y1Y=max(percaprev_pt(:,1)./percaprev_pt(:,2));
	set(ax(1),'YTick',round(([0:6]*(y1Y-y10)/6 +y10)*100)/100 );
	set(ax(1),'YLim', round(([0 6]*(y1Y-y10)/6 +y10)*100)/100 );
	set(ax(2),'YTick',round(([0:6]*(y2Y-y20)/6 +y20)*100)/100);
	set(ax(2),'YLim',round(([0 6]*(y2Y-y20)/6 +y20)*100)/100);
	set(gcf,'color','white');
	grid on;
	if (save_plots == 1) 
		saveas(gca,'rev_Pt_trans','eps2c'); 
		saveas(gca,'rev_Pt_trans.png'); 	
	end

	figure(6);
	h=plot([1:TT],percaprev_p0(:,1)./percaprev_p0(:,2) ,[1:TT], percaprev_p0(:,1));title('Relative Revenue Per Capita, P_0','FontSize',14);
	set(h,'LineWidth',2);legend('Location','North','P_0 y_a/y_m','P_0 y_a');
	set(gcf,'color','white');
	grid on;
	if (save_plots == 1) 
		saveas(gca,'rev_P0_trans','eps2c'); 
		saveas(gca,'rev_P0_trans.png'); 
	end


	figure(7);
	h=plot([1:TT],percaprev_pTT(:,1)./percaprev_pTT(:,2) ,[1:TT], percaprev_pTT(:,1));title('Relative Revenue Per Capita, P_T','FontSize',14);
	set(h,'LineWidth',2);legend('Location','North','P_T y_a/y_m','P_T y_a');
	set(gcf,'color','white');
	grid on;
	if (save_plots == 1) 
		saveas(gca,'rev_PTT_trans','eps2c'); 
		saveas(gca,'rev_PTT_trans.png'); 
	end

	figure(8);
	[ax,h1,h2]=plotyy([1:TT],price_path(:,2),[1:TT],Aa_path);title('Relative price of agricultural good','FontSize',14);
	set(h1,'LineWidth',2);set(h2,'LineWidth',2);legend('Location','North','P_t','A_a');
	set(gcf,'color','white');
	grid on;
	set(ax(1), 'YLim', [0.4 1.2]);
	set(ax(1), 'YTick', [0.4:0.2:1.2]);
	set(ax(2), 'YLim', [1.0 10.0]);
	set(ax(2), 'YTick', [1.0:3*0.75:10.0]);
	if (save_plots == 1) 
		saveas(gca,'agPrice','eps2c'); 
		saveas(gca,'agPrice.png');
	end

	figure(9);
	[ax,h1,h2]=plotyy([1:TT],1./price_path(:,2),[1:TT],Ym_path./Aa_path);title('APG  and Rel. Price (a la Alvarez-Cuadrado & Poschke)','FontSize',14);
	set(h1,'LineWidth',2);set(h2,'LineWidth',2,'color','r');legend('Location','North','1/P_t','y_m/A_a');
	set(gcf,'color','white');
	grid on;
	%set(ax(1), 'YLim', [0.8 1.6]);
	%set(ax(1), 'YTick', [0.8:0.2:1.6]);
	%set(ax(2), 'YLim', [0.6 1.4]);
	%set(ax(2), 'YTick', [0.6:0.2:1.4]);
	set(ax(2),'ycolor','r') ;
	if (save_plots == 1) 
		saveas(gca,'agPrice_APG','eps2c'); 
		saveas(gca,'agPrice_APG.png');
	end

	
	%%
	cd ~/Documents/CurrResearch/Devt/Computation

%end


%% put all of the countries on the same plot:
TTmax = max(c_nqtr);
ACPppath_ci = zeros(TTmax ,Ncountry);
ACPt_ci = zeros(TTmax ,Ncountry);
ACPindic_ci = zeros(TTmax ,Ncountry) == 0;
ACPtfpratio_ci = zeros(TTmax,Ncountry);
apg_ci = zeros(TTmax,Ncountry);

for ci=1:Ncountry
	cd(['trans_AaYm_results/calAa_' c_name{ci} 'Pa_linYm'])
	load(['trans_space_' c_name{ci} '.mat']);
	percaprev_pt = [price_path(:,2).*Aa_path'.*trans_path(:,1).^(mu-1) Ym_path'./(1-trans_path(:,1)-trans_path(:,2))];  % ag,man
	
	indic_data = c_time(:,ci)>0;
	yr0 = min(c_time(indic_data,ci)); nyr = sum(indic_data);
	time_qtr = interp1(c_time(indic_data,ci)-yr0,...
						c_time(indic_data,ci), 0:0.25:nyr-1);
	%c_time_qtr(1:c_nqtr(ci),ci) = time_qtr(1:end-1)';
	
	ACPindic_ci(:,ci) = c_time_qtr(:,ci)>0;
	
	ACPppath_ci(ACPindic_ci(:,ci),ci) = 1./price_path(:,2);
	ACPt_ci(ACPindic_ci(:,ci),ci) =  time_qtr(1:end-1)';
	ACPtfpratio_ci(ACPindic_ci(:,ci),ci)  = Ym_path./Aa_path;
	
	apg_ci(ACPindic_ci(:,ci),ci)  = percaprev_pt(:,2)./percaprev_pt(:,1);
	
	if strcmp(c_name{ci},'GER')
		apg_ci(1,ci) =  apg_ci(2,ci);
	end
	
	cd ~/Documents/CurrResearch/Devt/Computation
end

figure(1);
h=plot(	 ACPt_ci(ACPindic_ci(:,1),1), ACPppath_ci(ACPindic_ci(:,1),1),'r-.' , ACPt_ci(ACPindic_ci(:,1),1), ACPtfpratio_ci(ACPindic_ci(:,1),1),'r'  ...
	...	,ACPt_ci(ACPindic_ci(:,2),2), ACPppath_ci(ACPindic_ci(:,2),2),'g-.' , ACPt_ci(ACPindic_ci(:,2),2), ACPtfpratio_ci(ACPindic_ci(:,2),2),'g' ...
		,ACPt_ci(ACPindic_ci(:,3),3), ACPppath_ci(ACPindic_ci(:,3),3),'b-.' , ACPt_ci(ACPindic_ci(:,3),3), ACPtfpratio_ci(ACPindic_ci(:,3),3),'b' ...
		,ACPt_ci(ACPindic_ci(:,4),4), ACPppath_ci(ACPindic_ci(:,4),4),'k-.' , ACPt_ci(ACPindic_ci(:,4),4), ACPtfpratio_ci(ACPindic_ci(:,4),4),'k' ...
		);
set(h,'LineWidth',2);
set(gcf,'color','white');
grid on;
%legend('location','Northeast','Canada, P','Canada, TFP y_m/A_a','Germany, P','Germany, TFP y_m/A_a','UK, P','UK, TFP y_m/A_a','USA, P', 'USA,TFP y_m/A_a');
legend('location','Northeast','Canada, P','Canada, TFP y_m/A_a','UK, P','UK, TFP y_m/A_a','USA, P', 'USA,TFP y_m/A_a');
saveas(gca,'agPrice_reltfp_all','eps2c'); 
saveas(gca,'agPrice_reltfp_all.png');
	
	
figure(2)

h=plot(	 ACPt_ci(ACPindic_ci(:,1),1), ACPppath_ci(ACPindic_ci(:,1),1),'r-.' , ACPt_ci(ACPindic_ci(:,1),1), apg_ci(ACPindic_ci(:,1),1),'r'  ...
	...	,ACPt_ci(ACPindic_ci(:,2),2), ACPppath_ci(ACPindic_ci(:,2),2),'g-.' , ACPt_ci(ACPindic_ci(:,2),2), apg_ci(ACPindic_ci(:,2),2),'g' ...
		,ACPt_ci(ACPindic_ci(:,3),3), ACPppath_ci(ACPindic_ci(:,3),3),'b-.' , ACPt_ci(ACPindic_ci(:,3),3), apg_ci(ACPindic_ci(:,3),3),'b' ...
		,ACPt_ci(ACPindic_ci(:,4),4), ACPppath_ci(ACPindic_ci(:,4),4),'k-.' , ACPt_ci(ACPindic_ci(:,4),4), apg_ci(ACPindic_ci(:,4),4),'k' ...
		);
set(h,'LineWidth',2);
set(gcf,'color','white');
grid on;
%legend('location','Northeast','Canada, P','Canada, APG non-ag/ag','Germany, P','Germany, APG non-ag/ag','UK, P','UK, APG non-ag/ag','USA, P', 'USA, APG non-ag/ag');
legend('location','Northeast','Canada, P','Canada, APG non-ag/ag','UK, P','UK, APG non-ag/ag','USA, P', 'USA, APG non-ag/ag');
	saveas(gca,'agPrice_APG_all','eps2c'); 
	saveas(gca,'agPrice_APG_all.png');
