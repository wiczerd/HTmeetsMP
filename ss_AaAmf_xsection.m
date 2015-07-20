% solves transition paths from low to high steady states

% This program computes the transition for Harris-Todaro Meets
% Mortensen-Pissarides

% This program experiments with slower transitions.

cd ~/Documents/CurrResearch/Devt/Computation

global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau converged
global pfixed_devd Pa_devd
% load everything else (calibration and such):
% load trans_AaAmf/calAa_linYmAmf/trans_space_calAa_linYmAmf.mat
% load trans_AaAmf/trans_space_USA.mat 
load trans_AaYm_results/calAa_USAPa_linYm/trans_space_USA.mat

Amf_undevd = Amf_devd;
Pa_devd = wcPa_devd(2);

data = csvread('../Data/apg_u_na_gdp_wdi.csv',1,0);

Ncountry = size(data,1);

resid_ci =ones(Ncountry,3);
resid_AmfAa_ci=ones(Ncountry,3);
x_ci	=ones(Ncountry,3);
x_AmfAa_ci =ones(Ncountry,3);
x0_ci	=ones(Ncountry,3);

maxinitset = 20;
fval_x0_resid_x=zeros(maxinitset,10);


% cycle over pfixed
for pfixed_devd = 0:1

%% cycle over countries 
for ci =1:Ncountry

		Na_ci_target = data(ci,4);
		APG_ci_target = data(ci,2);
		u_ci_target = data(ci,3)/100;
		AmfAa_bnds = [0 0 0; 30 30 30];
		minresid = 100;minresid_i=1;
		for initset = 1:maxinitset
			if initset == 1
				x0 = [Amf_undevd, Aa_undevd Ym_undevd];
				Ym0 = Ym_undevd;
			elseif initset ==2
				x0 = [Amf_devd, Aa_devd  Ym_devd];
				Ym0 = Ym_devd;
			else
				x0 = rand([1 3]).*(AmfAa_bnds(2,:)-AmfAa_bnds(1,:)) + AmfAa_bnds(1,:);
			end
			%minimize APG gap
			cal_AmfAaYm_fn = @(AmfAaYm) calls_xsec_AmfAa(AmfAaYm,Na_ci_target,u_ci_target,APG_ci_target);
			[AmfAaYm1,fval_cal1,resid1,exitflag_cal1,out1] = lsqnonlin(cal_AmfAaYm_fn,x0,AmfAa_bnds(1,:),AmfAa_bnds(2,:),optsoff);
			fval_x0_resid_x(initset,:) = [fval_cal1,x0,resid1,AmfAaYm1];
			if(fval_cal1<minresid-1e-5 && exitflag_cal1 >= 0)
				minresid=fval_cal1;
				minresid_i = initset;
			elseif(initset - minresid_i >= 3 +floor(initset/3) )
				fval_cal1 = fval_x0_resid_x(minresid_i,1);
				x0 = fval_x0_resid_x(minresid_i,2:4);
				resid1 = fval_x0_resid_x(minresid_i,5:7);
				AmfAaYm1 = fval_x0_resid_x(minresid_i,8:10);

				break;
			end
		end
		
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
		[AmfAa1,fval_cal2,resid2,exitflag_cal1,out1] = lsqnonlin(cal_AmfAa_fn ,x_ci(1:2),AmfAa_bnds(1,:),AmfAa_bnds(2,:),optsoff);
		resid2 = cal_AmfAaYm_fn([AmfAa1 1]);
		resid_AmfAa_ci(ci,:) = resid2;
		x_AmfAa_ci(ci,:) =[AmfAa1 1];
		if exitflag_cal1<1
			disp(exitflag_cal1)
			disp(fval_cal1)
		end
		
	
end % loop over country
			%APG_ci_target - 
APG_ci	= data(:,2) - resid_ci(:,3)*100;
u_ci	= data(:,3)/100 - resid_ci(:,2);
Na_ci	= data(:,4) - resid_ci(:,1);
Pa_ci	= APG_ci./x_ci(:,2).*x_ci(:,3).*Na_ci.^(1-mu);


APG_AmfAa_ci	= data(:,2) - resid_AmfAa_ci(:,3)*100;
u_AmfAa_ci		= data(:,3)/100 - resid_AmfAa_ci(:,2);
Na_AmfAa_ci		= data(:,4)/100 - resid_AmfAa_ci(:,1);


cd xsec
save(['xsec_ss_cal_pfixed_devd' num2str(pfixed_devd) '.mat']);
cd ../
end



%% make some cross-sectional plots:

cd xsec
pfixed_devd = 1;
load(['xsec_ss_cal_pfixed_devd' num2str(pfixed_devd) '.mat']);


figure(1);
h=scatter(Pa_ci,Na_ci, 'filled','d');title('Labor in Agriculture, Calibrating to World Prices','FontSize',14);
xlabel('Relative Price of Agricultural Good','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['Na_ci_Pa_ci' num2str(pfixed_devd) ],'eps2c');
	saveas(gca,['Na_ci_Pa_ci' num2str(pfixed_devd) '.png']); 
end

figure(2);
h=scatter(Na_ci,APG_ci, 'filled','d');title('APG and Agricultural Labor, Calibrating to World Prices','FontSize',14);
xlabel('Fraction in agriculture','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['APG_ci_Na_ci' num2str(pfixed_devd) ],'eps2c');
	saveas(gca,['APG_ci_Na_ci' num2str(pfixed_devd) '.png']); 
end

figure(3);
h=scatter(u_ci,APG_ci, 'filled','d');title('APG and Unemployment, Calibrating to World Prices','FontSize',14);
xlabel('Unemployment Rate','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['APG_ci_u_ci' num2str(pfixed_devd)],'eps2c');
	saveas(gca,['APG_ci_u_ci' num2str(pfixed_devd) '.png']); 
end


figure(3);
h=scatter(x_ci(:,1),APG_ci, 'filled','d');title('APG and Matching Frictions, Calibrating to World Prices','FontSize',14);
xlabel('Matching function Efficiency','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['APG_ci_Amf_ci' num2str(pfixed_devd)],'eps2c');
	saveas(gca,['APG_ci_Amf_ci' num2str(pfixed_devd) '.png']); 
end

Amf_ci  = x_ci(:,1);
figure(4);
h=scatter(Amf_ci,APG_ci, 'filled','d');title('APG and Matching Efficiency, Calibrating to World Prices','FontSize',14);
xlabel('Matching Function Efficiency','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['APG_ci_Amf_ci' num2str(pfixed_devd)],'eps2c');
	saveas(gca,['APG_ci_Amf_ci' num2str(pfixed_devd) '.png']); 
end



pfixed_devd = 0;
load(['xsec_ss_cal_pfixed_devd' num2str(pfixed_devd) '.mat']);

figure(1);
h=scatter(Pa_ci,Na_ci, 'filled','d');title('Labor in Agriculture, Calibrating to Local Prices','FontSize',14);
xlabel('Relative Price of Agricultural Good','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['Na_ci_Pa_ci' num2str(pfixed_devd) ],'eps2c');
	saveas(gca,['Na_ci_Pa_ci' num2str(pfixed_devd) '.png']); 
end

figure(1);
h=scatter(Pa_ci,x_ci(:,2), 'filled','d');title('TFP in Agriculture, Calibrating to Local Prices','FontSize',14);
xlabel('Relative Price of Agricultural Good','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['Aa_ci_Pa_ci' num2str(pfixed_devd) ],'eps2c');
	saveas(gca,['Aa_ci_Pa_ci' num2str(pfixed_devd) '.png']); 
end


figure(1);
h=scatter(Pa_ci,x_ci(:,2)./x_ci(:,3), 'filled','d');title('TFP in Agriculture Relative to Non-Ag, Calibrating to Local Prices','FontSize',14);
xlabel('Relative Price of Agricultural Good','FontSize',14);
set(gcf,'color','white');
grid on;
if (save_plots == 1) 
	saveas(gca,['AaYm_ci_Pa_ci' num2str(pfixed_devd) ],'eps2c');
	saveas(gca,['AaYm_ci_Pa_ci' num2str(pfixed_devd) '.png']); 
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



cd ../
