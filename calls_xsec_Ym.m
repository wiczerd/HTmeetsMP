function APGresid = calls_xsec_Ym(Ym_here,AmfAa_0, AmfAa_bnds, Natarget,utarget,APGtarget)
global Ym
Ym = Ym_here;

optsoff = optimset('Display','off');
f_AmfAa = @(AmfAa) calls_xsec_AmfAa(AmfAa,Natarget,utarget);
[AmfAa1,fval_cal1,resid1,exitflag_cal1] = lsqnonlin(f_AmfAa,AmfAa_0,AmfAa_bnds(1,:),AmfAa_bnds(2,:),optsoff);
APGresid = calls_xsec_AmfAa(AmfAa1,Natarget,utarget,APGtarget);
% some error stuff if exit flag doesn't work
	if fval_cal1>1e-4 || exitflag_cal1<0
		disp(['did not solve interior' num2str(fval_cal1)]);
	end
end

