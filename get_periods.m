function [period,coi]=get_periods(frq,coi)
period=flip(1./frq);
coi=1./coi;
end

