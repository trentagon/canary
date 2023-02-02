%% Constants [no choice permitted]

in.amu = 1.660538782e-27; %kilograms in an amu
in.myr2sec = 31557600e+6; %seconds in 1 million years
in.mars_r = 3389.5e+3; %Mars radius in meters
in.mars_g = 3.71; %Mars surface gravity in meters sec^-2

%Present day composition of Mars atmosphere
%https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
%These 5 species make up 99.85% of the atmosphere
in.pres_co2 = 0.951;
in.pres_n2  = 0.0259;
in.pres_ar  = 0.0194;
in.pres_o2  = 0.0016;
in.pres_co  = 0.0006;

%% Free parameters

in.file_description = ['This input file is used to recreate figure 5 from the main manuscript.'];

%Start and end time
in.t_start   = 700; % [Myr] 700 Myr = 3.8 Ga (time starts at 4.5 Ga)
in.t_end     = 4500; % [Myr] 4500 Myr = 0 Ga = Today

%Hold a species at a constant pressure during evolution? 1 = yes, 0 = no.
in.static.c  = 0;
in.static.n  = 0;
in.static.ar = 0;

in.zt = 0.2132; %Distance from homopause to exobase divided by temperature, matters for diffusion
in.f_sputtering = 0.5022; %Sputtering multiplier, applied to all species

in.euv_uncertainty = 0; %Uncertainty in the euv flux evolution, [-0.1,0.1]
in.lym_uncertainty = 0; %Uncertainty in the lyman flux evolution, [-0.1,0.1]

in.euv_pl = 2.5307; %Power law index for N2 photochemical loss via euv flux
in.lym_pl = 2.9645; %Power law index for CO2 photochemical loss via lyman flux

in.f_outgassing = 1.5478; %Outgassing multiplier, applied to all species

in.f_photoC = 1.9524; %Photochemical multiplier for CARBON only
in.f_photoN = 0.7970; %Photochemical multiplier for NITROGEN only

in.f_ce = 0.3997; % Fraction of 40Ar released from crust

in.collapse = "obliquity"; %How to handle collapse of CO2 to the polar caps
%"obliquity": factor in obliquity variations
%"static": collapse under total pressure threshold
%"nothing": do not consider collapse

in.pc_col = 6; %[mbar] how much CO2 in the atmosphere during collapse?

in.comet = "no"; %Is there a comet?
in.comet_time = 3500; %When does it hit? [Myr]

%% Select chronology and load outgassing profile [All species]

in.chron_str = "hart"; %Which chronology to use? Hartmann 2005 or Ivanov 2001
%"hart" or "ivan"

chron_ivan = in.t_end - 1000*[3.86,3.74,3.65,3.46,1.45,0.387]; %Ivanov 2001
chron_hart = in.t_end - 1000*[3.85,3.57,3.40,3.00,0.880,0.235]; %Hartmann 2005

%Select chronology and corresponding crustal production profile
switch in.chron_str
    case "hart"
        chron = chron_hart;
        oo = load('og_hart.mat');
        in.og_profile = oo.og_profile;
    case "ivan"
        chron = chron_ivan;
        oo = load('og_ivan.mat');
        in.og_profile = oo.og_profile;
end

in.og_convert = 1000^3 * in.mars_g * 0.01 / (4 * pi * (in.mars_r)^2); %kg -> mbar
%note: og_profile is the crustal production rate in units of km^3/Myr
%inclusion of 1000^3 to convert og_profile to m^3/Myr

in.crust_density = 2900; %kg/m^3

in.og_profile(1,:) = in.f_outgassing * in.og_profile(1,:) * in.og_convert * in.crust_density;
%This currently gives the mbar/Myr of outgassing of all crustal species
%--> Multiply by a concentration [by mass] to get the outgassing profile
%[mbar/Myr] of a given species

%Static magma concentrations
%mag = concentration [by mass] of species in outgassing
in.mag_co2 = 9.654e-6; %Depends heavily on redox state of mantle
in.mag_n2 = 1.9E-6; %In overleaf
in.mag_ar36 = 3.46E-11; %In overleaf
in.mag_ar36to38 = 5.305; % ratio of ar36 to ar40 (see S&J 2016)
in.mag_ar38 = in.mag_ar36/in.mag_ar36to38; %Using ratio of ar36 to ar40 (see S&J 2016)

%Isotopic composition of mantle
in.dc_mantle  = 1 + (-0.025); %upper limit = -0.015, lower limit = -0.035
in.dn_mantle  = 1 + (-.03);   %upper limit = .3, lower limit = -0.03

%Starting delta values
in.dc_0  = in.dc_mantle;
in.dn_0  = in.dn_mantle;
in.ar36to38_start = in.mag_ar36to38;
in.ar40to36_start = 626;

%% Ar40 and K40 mantle calculations [Argon]
% Need to start from K40 concentration at 4.4 Ga and integrate to 
% get interior values and erosion rate starting point for model

mantle_vol = (4*pi/3)*(3340^3 - 1830^3)*10^9; % volume of mantle in m3
mantle_density = 3550; % kg/m3
in.mantle_mass = mantle_vol*mantle_density; % mass of mantle in kg

in.lambda_k40 = 5.543e-4; % total rate of K40 decay [Myr-1]
in.beta_ar40 = 0.1048; % branching ratio of K40 decay into Ar40
in.cx = 5; % how much K goes into the crust when outgassed?

switch in.chron_str
    case "hart"
        oo_ar = load('og_hart_from_44.mat');
        in.og_profile_ar = oo_ar.og_profile;
    case "ivan"
        oo_ar = load('og_ivan_from_44.mat');
        in.og_profile_ar = oo_ar.og_profile;
end

[in.ar40m_profile, in.crustal_erosion_rate] = calculate_40Ar_profiles(in);

%% Photochemical Loss [Carbon and Nitrogen]

%Escaping carbon particles per second
in.photo_flux_c = 4.68E+23; % Lo calculations
%in.photo_flux_c = 7.9E+23; % Hu et al. 2015

in.alph_photo_c = 0.7305; %Weighted average fractionation factor of carbon loss via all photochemical reactions
%From Lo calculations

%Escaping nitrogen particles per second: See Hu and Thomas 2022
in.photo_fluxPD_n = 1.4438e+23 * 0.41; %Photodissociation
in.photo_fluxDR_n = 2.8874e+23; %Dissociative Recombination
in.photo_fluxOP_n = 2.5411e+23; %Other processes

%Fractionation factors of nitrogen loss via photochemical reactions
in.alph_photoPD_n = 0.29; %Photodissociation
in.alph_photoDR_n = 0.58; %Dissociative Recombination
in.alph_photoOP_n = 1; %Other processes


%% Ion Loss [Carbon and Nitrogen]

in.f_ion = 1; %Ion multiplier
in.cphi0 = 8.66E+22; %Present day ion escape flux of CO2+ [particles/sec]
in.ion_ratio = 0.1;  %Ratio of N2+ to CO2+ at 160km, measured by MAVEN, Wither+15, Bougher+15


%% Carbonate and Nitrate Deposition [Carbon and Nitrogen]

%Carbonate Deposition-------
trans_time = 1000; %[Myr] Time of transition from early to late deposition rate
%Range = 1000 - 1500 Myr (= 3.5 - 3 Ga)

in.alph_carbonate = 1.01; %Fractionation factor of carbonate deposition
% Shallow subsurface aquifers: 1.01, and trans_time = 1500
% Open water systems: 1.06, and trans_time = 1000

early_total_c = 900; %[mbar] Amount of carbonate deposition before transition
late_total_c  = 7; %[mbar] Amount of carbonate deposition after transition (max = 7 mbar, Hu et al. 2015)

early_rate_c = early_total_c/abs(in.t_start - trans_time); %Convert total amount to rate
late_rate_c = late_total_c/abs(in.t_end - trans_time);

%Smooth the step function into a sigmoid
carb_x = in.t_start:0.01:in.t_end;
in.carbonate_profile(1,:) = early_rate_c + (sigmf(carb_x,[0.025, trans_time])*(late_rate_c - early_rate_c));
in.carbonate_profile(2,:) = carb_x;

%Nitrate Deposition-------
in.alph_nitrate = 1.01; %Fractionation factor of nitrate deposition

nit_massratio = 14/62; %Mass ratio in crust: N/NO3-
nit_conv = in.mars_g * 0.01; %kg(N)/ m^2 Myr -> mbar/Myr

%Amazonian values
nit_am_rho  = 2000; % [kg/m3] crust density
nit_am_dep  = 10; % [m] depth of crustal deposition
nit_am_conc = 180E-6; % [wt%] concentration of nitrate in crust

%Values in the Hesperian and Noachian (same units and meaning as Amazonian)
nit_nh_rho  = 3000;
nit_nh_dep  = 50; %d_nit: free parameter
nit_nh_conc = 300E-6;

early_total_n = nit_conv*nit_massratio*nit_nh_rho*nit_nh_dep*nit_nh_conc;
late_total_n = nit_conv*nit_massratio*nit_am_rho*nit_am_dep*nit_am_conc;

early_rate_n = early_total_n/abs(in.t_start - chron(4));%Convert total amount to rate
late_rate_n = late_total_n/abs(in.t_end - chron(4));

%Smooth the step function into a sigmoid
nit_x = in.t_start:0.01:in.t_end;
in.nitrate_profile(1,:) = early_rate_n + (sigmf(nit_x,[0.025, chron(4)])*(late_rate_n - early_rate_n));
in.nitrate_profile(2,:) = nit_x;

%% Load obliquity data if required

switch in.collapse
    case "obliquity"
        load("laskar_gaussian.mat"); %Load in obliquity probability distributions
        in.obliquity_data.profile = laskar_gaussian;
        
        load("forget_criticals.mat"); %Load in critical obliquity of collapse as a function of ptot
        in.obliquity_data.c_i = forget_profile;
    case "static"
        in.static_collapse_threshold = 500; %[mbar] total pressure threshold of CO2 collapse
    case "nothing"
        %Do nothing
end

%% Comet

in.par_comet = 0; % [mbar] Reset Ar to this pressure
in.dar_comet = 1 + (0); %Reset 1 (solar value is reference) (Look at Balsiger et al. 2015?)
in.impacted = 1; % [Do not change] Switch that turns to 0 when the comet impacts

%% Interplanetary Dust Particles [Argon]

in.f_idp = 1; %Multiplier of interplanetary dust particles
in.idp_rate_tot = in.f_idp*2.878e-10; %[mbar] Rate of deposition (Flynn 1997) 

idp_36to38 = 5; %36Ar/38Ar = 5
in.idp_rate_36 = in.idp_rate_tot/(1+(1/idp_36to38));
in.idp_rate_38 = in.idp_rate_tot/(1+(idp_36to38));

%% Technical specifications

%Timestepping
in.ptol_fast = 0.15; %Maximum percent change to pressure via fastest source or sink
in.ptol_total = 0.01; %Maximum percent change to pressure via all processes combined
in.k_max = 1; %Maximum allowable timestep
in.k_min = 0.1; %Minimum allowable timestep

%Negative atmosphere handling
in.min_pressure = 0.0001; %[mbar] Set the pressure to this value if it is going to be negative
in.min_pressure_ar = 1e-12; %[mbar] special lower pressure for the super low ar36 and ar38
%pCO2 cannot be zero because sputtering divides by it for scaling

%Pre prescribe an evolution?
in.prescribe_c = 0;
in.preData_c = 0;

in.prescribe_n = 0;
in.preData_n = 0;

in.prescirbe_ar = 0;
in.preData_ar = 0;





