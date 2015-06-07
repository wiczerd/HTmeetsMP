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


data = csvread('../Data/apg_u_na_gdp_wdi.csv',1,0);

converged = 0;

%% transition backwards 
%  First, transition using the approximation that unemployment is the steady state.  
%  This doesn't matter explicitly for any of the choices, just market clearing.

Aa = Aa_devd;

%%  set up new sequences

%set up the stunted growth paths in a range, and then put them together for
%the x-section

Ncountry = size(data,1);

resid_ci =ones(Ncountry,3);
resid_AmfAa_ci=ones(Ncountry,3);
x_ci	=ones(Ncountry,3);
x_AmfAa_ci =ones(Ncountry,3);
x0_ci    =ones(Ncountry,3);
abar_0 = abar;
for ci =1:Ncountry

		Na_ci_target = data(ci,4)/100;
		APG_ci_target = data(ci,2);
		u_ci_target = data(ci,3)/100;
		if Na_ci_target >.5
			x0 = [Amf_undevd, Aa_undevd Ym_undevd];
			Ym0 = Ym_undevd;
		else
			x0 = [Amf_devd, Aa_devd  Ym_devd];
			Ym0 = Ym_devd;
		end
		%x0(2) = APG_ci_target*x0(1);
		
		%minimize APG gap
		AmfAa_bnds = [0 0 0; 30 30 30];
		cal_AmfAaYm_fn = @(AmfAaYm) calls_xsec_AmfAa(AmfAaYm,Na_ci_target,u_ci_target,APG_ci_target);
		[AmfAaYm1,fval_cal1,resid1,exitflag_cal1,out1] = lsqnonlin(cal_AmfAaYm_fn,x0,AmfAa_bnds(1,:),AmfAa_bnds(2,:),optsoff);
		
		resid_ci(ci,:) = resid1 ;
		x_ci(ci,:) =AmfAaYm1;
		x0_ci(ci,:) = x0 ;
		if exitflag_cal1<1
			disp(exitflag_cal1)
			disp(fval_cal1)
		end
		%don't touch APG
		AmfAa_bnds = AmfAa_bnds(:,1:2);
		cal_AmfAa_fn = @(AmfAa) calls_xsec_AmfAa(AmfAa,Na_ci_target,u_ci_target);
		[AmfAa1,fval_cal1,resid1,exitflag_cal1,out1] = lsqnonlin(cal_AmfAa_fn ,AmfAaYm1(1:2),AmfAa_bnds(1,:),AmfAa_bnds(2,:),optsoff);
		resid1 = cal_AmfAaYm_fn([AmfAa1 1]);
		resid_AmfAa_ci(ci,:) = resid1 ;
		x_AmfAa_ci(ci,:) =[AmfAa1 1];
		if exitflag_cal1<1
			disp(exitflag_cal1)
			disp(fval_cal1)
		end
		
	
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