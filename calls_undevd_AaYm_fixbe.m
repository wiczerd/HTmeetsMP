function calresid_v = calls_undevd_AaYm_fixbe(Aabeabar,Natarget,Patarget)
% the objective for a calibration to target Na, uss and solve for prices

global cbar abar Aa beta eta Ym lambda kappa theta Amf mu alpha be tau


Aa  = (Aabeabar(1));
%be = (Aabeabar(2));
abar= (Aabeabar(2));
options = optimset('Display','off');

tau=0.0;

pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);

[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(.5*pi/Ym-pi/2) log(.5)],options);
[excess,theeconomy] = sol_wcPa_ss([(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))]);

if(exitflag <0)
	calresid_v = 10*ones(3,1);
else
	% loop on tau
	tauH = .1;tauL=0.;
	for itertau = 1:100
		tau = 0.5*tauH+0.5*tauL;
		[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,logssp,options);
		wcPa_ss =[(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
		% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
		[excess,theeconomy] = sol_wcPa_ss([wcPa_ss]);
		budget_def = be*theeconomy(2) - wcPa_ss(1)*tau*(1-theeconomy(2)-theeconomy(1));
		if(abs(budget_def)<1e-6 || (tauH-tauL)<1e-6)
			break;
		elseif (budget_def < 0)
			tauH=tau;
		elseif(budget_def > 0)
			tauL=tau;
		end
	end

	calresid_v(1)	= Natarget - theeconomy(1);
%	calresid_v(2)	= utarget - theeconomy(2)/(1-theeconomy(1));
	calresid_v(2)	= Patarget - wcPa_ss(2);
	% fall off a cliff for unemployment if everyone is in agriculture
	if(theeconomy(1)>=1) calresid_v(2) = 50; end%.05*(1/Amf); end
end	
