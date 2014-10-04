function calresid_v = calls_undevd_beabar(beabar,Pa,w0,Natarget,utarget,trans_path)
% the objective for a calibration to target Na, uss and solve for prices

global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau

be = (beabar(1));
abar= (beabar(2));
%abar = beabar;
options = optimset('Display','off');

tau=0.0;

pos_solwc = @(wc) sol_wc((atan(wc)+pi/2)*Ym/pi,Pa,trans_path(2,:),trans_path(1,2)); % should I use trans_path(1,2) instead of utarget*(1-Natarget)
[tanw, fval,exitflag,output,J] = fsolve(pos_solwc,tan(w0*pi/Ym-pi/2) ,options);
[excess,theeconomy] = sol_wc((atan(tanw)+pi/2)*Ym/pi,Pa,trans_path(2,:),trans_path(1,2));
%pos_solwcAa = @(wcAa) sol_wcAa([(atan(wcAa(1))+pi/2)*Ym/pi exp(wcAa(2))],Pa,trans_path(2,:),trans_path(1,2)); % should I use trans_path(1,2) instead of utarget*(1-Natarget)
%[tanw, fval,exitflag,output,J] = fsolve(pos_solwcAa,[tan(w0*pi/Ym-pi/2) log(Aa)],options);

%[excess,theeconomy] = sol_wcAa([(atan(tanw(1))+pi/2)*Ym/pi exp(tanw(2))],Pa,trans_path(2,:),trans_path(1,2));

if(exitflag <0)
	calresid_v = 10;
else
	% loop on tau
	tauH = .1;tauL=0.;
	for itertau = 1:100
		tau = 0.5*tauH+0.5*tauL;
		[tanw, fval,exitflag,output,J] = fsolve(pos_solwc,tan(w0*pi/Ym-pi/2),options);
		[excess,theeconomy] = sol_wc((atan(tanw)+pi/2)*Ym/pi,Pa,trans_path(2,:),trans_path(1,2));
		wc = (atan(tanw)+pi/2)*Ym/pi ;
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

	calresid_v(1)	= Natarget - theeconomy(1);
	calresid_v(2)	= utarget  - theeconomy(2)/(1-theeconomy(1));

	% fall off a cliff for unemployment if everyone is in agriculture
%	if(theeconomy(1)>=1) calresid_v(2) = 50; end%.05*(1/Amf); end
end	
