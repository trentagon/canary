function [pc_new, pn_new, par36_new, par38_new, par40_new, flag] = advancePressure(sv,in,rates)

%Forward Euler advance:
%First, check if the species is held static
if in.static.c
    pc_new = sv.pc;
else
    pc_new = sv.pc + sum(rates.c)*sv.k;
end

if in.static.n
    pn_new = sv.pn ;
else
    pn_new = sv.pn + sum(rates.n)*sv.k;
end

if in.static.ar
    par36_new = sv.par36;
    par38_new = sv.par38;
    par40_new = sv.par40;
else
    par36_new = sv.par36 + sum(rates.ar36)*sv.k;
    par38_new = sv.par38 + sum(rates.ar38)*sv.k;
    par40_new = sv.par40 + sum(rates.ar40)*sv.k;
end

%Handle negative atmospheres
%{
A negative pressure can occur because nitrate and carbonate deposition are
fixed rates, and do not depend on the mixing ratios. All other loss
processes will adjust 
%}

%Set negative flags to 0
flag.c = 0;
flag.n = 0;
flag.ar36 = 0;
flag.ar38 = 0;
flag.ar40 = 0;

%Fix the pressure
if pc_new < in.min_pressure
    pc_new = in.min_pressure;
    flag.c = 1;
end

if pn_new < in.min_pressure
    pn_new = in.min_pressure;
    flag.n = 1;
end

if par36_new < in.min_pressure_ar
    par36_new = in.min_pressure;
    flag.ar36 = 1;
end

if par38_new < in.min_pressure_ar
    par38_new = in.min_pressure;
    flag.ar38 = 1;
end

if par40_new < in.min_pressure
    par40_new = in.min_pressure;
    flag.ar40 = 1;
end

end