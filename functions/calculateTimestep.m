function k = calculateTimestep(sv,in,rates)
%{
This function calculates the maximum allowable timestep for the model at a 
given time. The principle behind this calculation is that at each timestep,
we require that the pressure of the fastest changing species only change by 
the percentage 'in.ptol_total'. We're using Forward Euler, so this can 
be done exactly.
%}

%Find the fastest rate for each species
fastest_c = max(abs(rates.c));
fastest_n = max(abs(rates.n));
fastest_ar36 = max(abs(rates.ar36));
fastest_ar38 = max(abs(rates.ar38));
fastest_ar40 = max(abs(rates.ar40));

%Find max timestep based on change threshold and fastest rate
k_fastest = in.ptol_fast*min([sv.pc/fastest_c, sv.pn/fastest_n, sv.par36/fastest_ar36, sv.par38/fastest_ar38, sv.par40/fastest_ar40]);

% k_fastest = in.ptol_fast*min([sv.pc/fastest_c, sv.pn/fastest_n, sv.par40/fastest_ar40]);


%Find the total rate for each species
total_c = abs(sum(rates.c));
total_n = abs(sum(rates.n));
total_ar36 = abs(sum(rates.ar36));
total_ar38 = abs(sum(rates.ar38));
total_ar40 = abs(sum(rates.ar40));

%Find max timestep based on change threshold and total rate
k_total = in.ptol_total*min([sv.pc/total_c, sv.pn/total_n, sv.par36/total_ar36, sv.par38/total_ar38, sv.par40/total_ar40]);

% k_total = in.ptol_total*min([sv.pc/total_c, sv.pn/total_n, sv.par40/total_ar40]);


%Take minimum of the two possible timesteps and the maximum allowable step
k = min([k_fastest, k_total, in.k_max]);

%Take maximum of k and the minimum allowable step
k = max([k, in.k_min]);

end