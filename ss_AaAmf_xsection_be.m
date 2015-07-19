% solves transition paths from low to high steady states

% This program computes the transition for Harris-Todaro Meets
% Mortensen-Pissarides

% This program experiments with slower transitions.

cd ~/Documents/CurrResearch/Devt/Computation

global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau converged
global pfixed_devd Pa_devd
% load everything else (calibration and such):
load trans_AaAmf/calAa_linYmAmf/trans_space_calAa_linYmAmf.mat
%load trans_AaAmf/trans_space_USA.mat 

Pa_devd = wcPa_devd(2);

pfixed_devd  = 0;

data = csvread('../Data/apg_u_na_gdp_wdi.csv',1,0);

Ncountry = size(data,1);

resid_ci =ones(Ncountry,3);
resid_beAa_ci=ones(Ncountry,3);
x_ci	=ones(Ncountry,3);
x_beAa_ci =ones(Ncountry,3);
x0_ci	=ones(Ncountry,3);


%% cycle over countries 

for ci =1:Ncountry

		Na_ci_target = data(ci,4);
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
		beAa_bnds = [0 0 0; Ym 30 30];
		cal_beAaYm_fn = @(beAaYm) calls_xsec_beAa(beAaYm,Na_ci_target,u_ci_target,APG_ci_target);
		[beAaYm1,fval_cal1,resid1,exitflag_cal1,out1] = lsqnonlin(cal_beAaYm_fn,x0,beAa_bnds(1,:),beAa_bnds(2,:),optsoff);
		
		resid_ci(ci,:) = resid1 ;
		x_ci(ci,:) =beAaYm1;
		x0_ci(ci,:) = x0 ;
		if exitflag_cal1<1
			disp(exitflag_cal1)
			disp(fval_cal1)
		end
		%don't touch APG
		beAa_bnds = beAa_bnds(:,1:2);
		cal_beAa_fn = @(beAa) calls_xsec_beAa(beAa,Na_ci_target,u_ci_target);
		[beAa1,fval_cal2,resid2,exitflag_cal1,out1] = lsqnonlin(cal_beAa_fn ,x0(1:2),beAa_bnds(1,:),beAa_bnds(2,:),optsoff);
		resid2 = cal_beAaYm_fn([beAa1 1]);
		resid_beAa_ci(ci,:) = resid2;
		x_beAa_ci(ci,:) =[beAa1 1];
		if exitflag_cal1<1
			disp(exitflag_cal1)
			disp(fval_cal1)
		end
	
end % loop over country
			%APG_ci_target - 
APG_ci	= data(:,2) - resid_ci(:,3)*100;
u_ci	= data(:,3)/100 - resid_ci(:,2);
Na_ci	= data(:,4) - resid_ci(:,1);
Pa_ci	= 1./APG_ci./x_ci(:,2).*x_ci(:,3).*Na_ci.^(1-mu);


APG_AmfAa_ci	= data(:,2) - resid_beAa_ci(:,3)*100;
u_AmfAa_ci		= data(:,3)/100 - resid_beAa_ci(:,2);
Na_AmfAa_ci		= data(:,4)/100 - resid_beAa_ci(:,1);


cd xsec
save(['xsec_ss_cal_be_pfixed_devd' num2str(pfixed_devd) '.mat']);
cd ../



%% make some cross-sectional plots:

cd xsec

pfixed_devd = 0;
load(['xsec_ss_cal_be_pfixed_devd' num2str(pfixed_devd) '.mat']);
mkdir('xsec_be')
cd xsec_be

figure(1);
h=scatter(Pa_ci,Na_ci, 'filled','d');title('Labor in Agriculture, Calibrating to Local Prices','FontSize',14);
xlabel('Relative Price of Agricultural Good','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['Na_ci_Pa_ci' num2str(pfixed_devd) ],'eps2c');
	saveas(gca,['Na_ci_Pa_ci' num2str(pfixed_devd) '.png']); 
end

figure(2);
h=scatter(Na_ci,APG_ci, 'filled','d');title('APG and Agricultural Labor, Calibrating to Local Prices','FontSize',14);
xlabel('Fraction in agriculture','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['APG_ci_Na_ci' num2str(pfixed_devd) ],'eps2c');
	saveas(gca,['APG_ci_Na_ci' num2str(pfixed_devd) '.png']); 
end

figure(3);
h=scatter(u_ci,APG_ci, 'filled','d');title('APG and Unemployment, Calibrating to Local Prices','FontSize',14);
xlabel('Unemployment Rate','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['APG_ci_u_ci' num2str(pfixed_devd)],'eps2c');
	saveas(gca,['APG_ci_u_ci' num2str(pfixed_devd) '.png']); 
end

Amf_ci  = x_ci(:,1);
figure(4);
h=scatter(Amf_ci,APG_ci, 'filled','d');title('APG and Matching Efficiency, Calibrating to Local Prices','FontSize',14);
xlabel('Matching Function Efficiency','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['APG_ci_Amf_ci' num2str(pfixed_devd)],'eps2c');
	saveas(gca,['APG_ci_Amf_ci' num2str(pfixed_devd) '.png']); 
end

figure(5);
h=scatter(data(:,3)/100, u_ci, 'filled','d');title('Unemployment Rate Fit, Calibrating to Local Prices','FontSize',14);
xlabel('Unemployment rate, data','FontSize',14);
ylabel('Unemployment rate, model','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['u_ci_cal_dat' num2str(pfixed_devd)],'eps2c');
	saveas(gca,['u_ci_cal_dat' num2str(pfixed_devd) '.png']); 
end


figure(6);
h=scatter(data(:,4),Na_ci, 'filled','d');title('Agricultural Labor Force Fit, Calibrating to Local Prices','FontSize',14);
xlabel('% Agriculture, data','FontSize',14);
ylabel('% Agriculture, model','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['Na_ci_cal_dat' num2str(pfixed_devd)],'eps2c');
	saveas(gca,['Na_ci_cal_dat' num2str(pfixed_devd) '.png']); 
end

figure(7);
h=scatter(data(:,2),APG_ci, 'filled','d');title('APG fit, Calibrating to Local Prices','FontSize',14);
xlabel('APG, data','FontSize',14);
ylabel('APG, model','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['APG_ci_cal_dat' num2str(pfixed_devd)],'eps2c');
	saveas(gca,['APG_ci_cal_dat' num2str(pfixed_devd) '.png']); 
end



cd ../../
