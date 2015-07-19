function [calresid_v] = calls_xsec_beAa(beAaYm,Natarget,utarget,APGtarget)
% the objective for a calibration to target Na, uss and solve for prices

global Aa Ym mu be tau 
global pfixed_devd Pa_devd

%Aa = (beAaYm(1));
Aa = (beAaYm(2));
be = (beAaYm(1));
if nargin>3
	Ym = (beAaYm(3));
else
	Ym = 1;
end
options = optimset('Display','off');

tau=0.0;
%pos_solwcPa = @(logwcPa) sol_wcPa_ss(exp(logwcPa));

pos_solwcPa = @(wcPa) sol_wcPa_ss([(atan(wcPa(1))+pi/2)*Ym/pi exp(wcPa(2))]);
wcPa_0 = [Ym/2. 1.];
[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,[tan(wcPa_0(1)*pi/Ym-pi/2) log(wcPa_0(2))],options);
[excess,theeconomy] = sol_wcPa_ss([(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))]);

if(exitflag <0)
%	[logssp, fval_r,resid_r,exitflag_r,ouput_r,J_r] = lsqnonlin(sol_wcPa_ss,[0.5 1.], [0. 0.], [Ym inf],options);
%	[excess,theeconomy] = sol_wcPa_ss([(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))]);
	calresid_v = 10.*ones(length(beAaYm),1);
	
else
	% loop on tau
	tauH = .1;tauL=0.;
	for itertau = 1:100
		tau = 0.5*tauH+0.5*tauL;
		[logssp, fval,exitflag,output,J] = fsolve(pos_solwcPa,logssp,options);
		wcPa_ss =[(atan(logssp(1))+pi/2)*Ym/pi exp(logssp(2))];
		% theeconomy{:} = {N_a, u, Q, J, Ve, Vu}
		[excess,theeconomy] = sol_wcPa_ss([wcPa_ss]);
		%absolute level of UI replacement
		%budget_def = be*theeconomy(2) - wcPa_ss(1)*tau*(1-theeconomy(2)-theeconomy(1));
		%propotional UI replacement
		budget_def = be*theeconomy(2) - tau*(1-theeconomy(2)-theeconomy(1));
		if(abs(budget_def)<1e-6 || (tauH-tauL)<1e-6)
			break;
		elseif (budget_def < 0)
			tauH=tau;
		elseif(budget_def > 0)
			tauL=tau;
		end
	end
	uhere = theeconomy(2)/(1-theeconomy(1));
	%percaprev_pt = [wcPa_ss(2).*Aa*theeconomy(1).^(mu-1) Ym/(1-uhere )];  % ag,man
	percaprev_pt = [wcPa_ss(2).*Aa*theeconomy(1).^(mu-1) Ym];  % ag,man
	%percaprev_pt = [Aa*theeconomy(1).^(mu-1) Ym/(1-theeconomy(2) )];  % ag,man
	if pfixed_devd  ==1
		APG = percaprev_pt(2)/percaprev_pt(1)/wcPa_ss(2)*Pa_devd;
	else
		APG = percaprev_pt(2)/percaprev_pt(1);
	end
	calresid_v(1)	= Natarget - theeconomy(1);
	calresid_v(2)	= utarget - uhere;
	if nargin>3
		APGresid		= APGtarget - APG;
		calresid_v(3) = APGresid/100;
	end
	if(theeconomy(1)>=1) calresid_v(2) = 50; end%.05*(1/Amf); end
	if(any(~isfinite(calresid_v)))
		if(~isfinite(calresid_v(1))) calresid_v(1) = Natarget - 1. ; end;
		if(~isfinite(calresid_v(2))) calresid_v(2) = utarget - 1. ; end;
		if nargin>3 && ~isfinite(calresid_v(3))
				calresid_v(3) = APGtarget - 20;
		end;
	end
end	

%disp(calresid_v)