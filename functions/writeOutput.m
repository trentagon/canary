function out = writeOutput(sv,rates,out)

out.pressures.c(end+1,1)  = sv.pc;
out.pressures.n(end+1,1)  = sv.pn;
out.pressures.ar36(end+1,1)  = sv.par36;
out.pressures.ar38(end+1,1)  = sv.par38;
out.pressures.ar40(end+1,1)  = sv.par40;

out.deltas.c(end+1,1)  = sv.dc;
out.deltas.n(end+1,1)  = sv.dn;
out.deltas.ar38(end+1,1)  = sv.dar38;
out.deltas.ar40(end+1,1)  = sv.dar40;

out.rates.c(end+1,:) = rates.c';
out.rates.n(end+1,:) = rates.n';
out.rates.ar36(end+1,:) = rates.ar36';
out.rates.ar38(end+1,:) = rates.ar38';
out.rates.ar40(end+1,:) = rates.ar40';

out.times.t(end+1,1) = sv.t;
out.times.k(end+1,1) = sv.k;

out.negative_flag.c(end+1,1)  = sv.negative_flag.c;
out.negative_flag.n(end+1,1)  = sv.negative_flag.n;
out.negative_flag.ar36(end+1,1) = sv.negative_flag.ar36;
out.negative_flag.ar38(end+1,1) = sv.negative_flag.ar38;
out.negative_flag.ar40(end+1,1) = sv.negative_flag.ar40;

out.collapse(end+1,1) = sv.fc;
end 