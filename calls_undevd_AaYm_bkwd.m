function calresid_v = calls_undevd_AaYm_bkwd(Aaabar,Pa,w0,Natarget,trans_path)
% the objective for a calibration to target Na, uss and solve for prices

global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau

Aa = Aaabar(1);
abar= Aaabar(2);
%abar = beabar;
options = optimset('Display','off');

tau=0.0;

pos_solwcPa = @(wcPa) sol_wcPa([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))],trans_path(2,:),trans_path(1,2)); % should I use trans_path(1,2) instead of utarget*(1-Natarget)
[tanw, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(w0*pi/Ym-pi/2) log(Pa)],options);

if(exitflag <0)
	calresid_v(1) = 10;
	calresid_v(2) = 10;
else
	% loop on tau
	tauH = .1;tauL=0.;
	for itertau = 1:100
		tau = 0.5*tauH+0.5*tauL;

		[tanw, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(w0*pi/Ym-pi/2) tanw(2)],options);
		[excess,theeconomy] = sol_wcPa([(atan(tanw(1))+pi/2)*Ym/pi exp(tanw(2))],trans_path(2,:),trans_path(1,2));
		wc = (atan(tanw(1))+pi/2)*Ym/pi ;
		
		% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
		budget_def = be*theeconomy(2) - wc*tau*(1-theeconomy(2)-theeconomy(1));
		if(abs(budget_def)<1e-6 || (tauH-tauL)<1e-6)
			break;
		elseif (budget_def < 0)
			tauH=tau;
		elseif(budget_def > 0)
			tauL=tau;
		end
	end

	calresid_v(2)	= Natarget - theeconomy(1);
	calresid_v(1)	= Pa - exp(tanw(2));

% 	
% 	abar = abar-5e-4;
% 	[tanw_m, fval_m,exitflag_m,output_m,J_m] = fsolve(pos_solwcPa,[tan(w0*pi/Ym-pi/2) tanw(2)],options);
% 	[excess_m,theeconomy_m] = sol_wcPa([(atan(tanw(1))+pi/2)*Ym/pi exp(tanw(2))],Pa,trans_path(2,:),trans_path(1,2));
% 	abar = abar+1e-3;
% 	[tanw_p, fval_p,exitflag_p,output_p,J_p] = fsolve(pos_solwcPa,[tan(w0*pi/Ym-pi/2) tanw(2)],options);
% 	[excess_p,theeconomy_p] = sol_wcPa([(atan(tanw(1))+pi/2)*Ym/pi exp(tanw(2))],Pa,trans_path(2,:),trans_path(1,2));
% 	abar = abar-5e-4;
% 	
% 	Jhere = - (theeconomy_p(1)-theeconomy_m(1))/1e-3;
	% fall off a cliff for unemployment if everyone is in agriculture
%	if(theeconomy(1)>=1) calresid_v(2) = 50; end%.05*(1/Amf); end
end	
