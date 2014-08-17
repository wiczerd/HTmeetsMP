% solves transition paths from low to high steady states

% This program computes the transition for Harris-Todaro Meets
% Mortensen-Pissarides

cd ~/Documents/CurrResearch/Devt/Computation

global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau

TT = 100;
save_plots =0;

cbar	= -0.6; % note this is the inverse because I changed the util function
abar	= 0.09; % this gets changed below in the calibration
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

%% solve for the steady state with the original parameters

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
% unemployment by changing abar and Amf

% to change the calibration target values change Na_target, u_target
Na_target = 0.05;
u_target  = 0.05;
cal_devd55 = @(abarAmf) cal_devd(abarAmf,Na_target,u_target);
[x,fval_cal,exitflag_cal,out] = fminsearch(cal_devd55,log([abar,Amf]));

pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);
[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)]);
wcPa_ss = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
[excess_devd,devd_economy] = sol_wcPa_ss(wcPa_ss);
devd_logssp = logssp;
Aa_devd = Aa;
Amf_devd = Amf;
Ym_devd  = Ym;

%%
% now calibrate it so that I get 90% in agriculture and 10% unemployment by
% manipulating Aa and Amf.  I will fix abar.


% to change the calibration target values change Na_target, u_target
Na_target = 0.90;
u_target  = 0.10;
cal_undevd9010 = @(AaYm) cal_undevd_AaYm(AaYm,Na_target,u_target);
[x,fval_cal,exitflag_cal,out] = fminsearch(cal_undevd9010,log([.1,Ym]));


pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);
[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)]);
wcPa_ss = [(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
[excess_undevd,undevd_economy] = sol_wcPa_ss(wcPa_ss);

undevd_logssp = logssp;
Aa_undevd = Aa;
Ym_undevd = Ym;
Amf_undevd= Amf;

%% transition backwards 
%  First, transition using the approximation that unemployment is the steady state.  
%  This doesn't matter explicitly for any of the choices, just market clearing.

Aa = Aa_devd;

%Here is where you can change the rate of change of the two paths, making
%rate change power over 1 makes for a slow transition and less than 1 is a
%fast transition

rt_chng_pwr	= 1;
Aa_path_chng	= linspace(0,1,TT).^rt_chng_pwr;
Aa_path		= (1-Aa_path_chng)*Aa_undevd+Aa_path_chng*Aa_devd;

rt_chng_pwr	= 5;
Ym_path_chng	= linspace(0,1,TT).^rt_chng_pwr;
Ym_path		= (1-Ym_path_chng)*Ym_undevd+Ym_path_chng*Ym_devd;

trans_path  = zeros(TT,size(devd_economy,2));
excess_path = zeros(TT,2);
price_path  = zeros(TT,2);


p0_trans = [linspace(exp(devd_logssp(2)),exp(undevd_logssp(2)),100);...
	   linspace((atan(devd_logssp(1))+pi/2)*Ym/pi,(atan(undevd_logssp(1))+pi/2)*Ym/pi,100)]';

trans_economy	= devd_economy;
trans_path(TT,:)= trans_economy;

trans_economy	= devd_economy;
trans_path(TT,:)= trans_economy;


% use steady state unemployment first
trans_path(:,2) = -1.0;

price_path(TT,:)= p0_trans(1,:);

for trans_iter =1:10

for t = TT-1:-1:1
	Aa = Aa_path(t);
	Ym = Ym_path(t);
	pos_solwcPa = @(wcPa) sol_wcPa([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))],trans_path(t+1,:),trans_path(t,2));
	tauH = 0.05; tauL=0.001;
	for itertau = 1:100
		tau = 0.5*tauH+0.5*tauL;
		%p0 = [tan(p0_trans(t,1)*pi/Ym-pi/2) log(p0_trans(t,2))];
		if(t>TT/2)
			p0 = devd_logssp;
		else
			p0 = undevd_logssp;
		end
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
		if(t>TT/2)
			p0 = devd_logssp;
		else
			p0 = undevd_logssp;
		end
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
	price_path(t,:) = wcPa_t;
end
trans_path(TT,:) = devd_economy;
price_dif = abs(price_path - price_path_back);
if(max(max(price_dif))<1e-6)
	break;
end

end


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

cd ~/Documents/CurrResearch/Devt/Computation