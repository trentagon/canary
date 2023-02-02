function [c_rates, n_rates, ar36_rates, ar38_rates, ar40_rates] = calculateRates(pCO2,sv,in)
%% Calculate the rates of the sources and sinks for all 3 species
%{ 
Inputs: pCO2, statevector, inputs
Outputs: 1 vector per species containing the rates of change due to sources
and sinks [mbar/myr]
%}


%% Establish frequently used variables

%Pressures
pc = pCO2; %pCO2 is its own input variable to make handling of collapse easier
pn = sv.pn;
par36 = sv.par36;
par38 = sv.par38;
par40 = sv.par40;

%Mixing Ratios
mr_bot = (pn/28)+(pc/44)+(par36/36)+(par38/38)+(par40/40);
mrc    = (pc/44)/mr_bot;
mrn    = (pn/28)/mr_bot;
mrar36 = (par36/36)/mr_bot;
mrar38 = (par38/38)/mr_bot;
mrar40 = (par40/40)/mr_bot;

%Conversions: particles per second --> mbar/myr
co2_conv  = 44*(in.myr2sec*in.amu*in.mars_g*0.01)/(4 * pi * in.mars_r^2); %this is CO2!! not C
n2_conv  = 28*(in.myr2sec*in.amu*in.mars_g*0.01)/(4 * pi * in.mars_r^2); %this is N2!! not N
ar36_conv = 36*(in.myr2sec*in.amu*in.mars_g*0.01)/(4 * pi * in.mars_r^2);
ar38_conv = 38*(in.myr2sec*in.amu*in.mars_g*0.01)/(4 * pi * in.mars_r^2);
ar40_conv = 40*(in.myr2sec*in.amu*in.mars_g*0.01)/(4 * pi * in.mars_r^2);
% ^^^ note the 0.01 is there to go from pascal to mbar

%Get the index in outgassing and deposition profiles that corresponds to the model time 
idx  = floor(100*(round(sv.t,2) - 700) + 1); 

%% Sputtering [CNAR]

euv_flux = power(sv.t/4500,-1.23 + in.euv_uncertainty); %Evolution of the EUV flux
sput_co2_pps = exp(-0.462*power(log(euv_flux),2)+5.0866*log(euv_flux)+53.49); %sputtering of CO2 in particles per second

%Dilution factor
dilute = 1 + (mrn/mrc)*exp(.446*16*in.zt) + (mrar36/mrc)*exp(.446*8*in.zt) + (mrar38/mrc)*exp(.446*6*in.zt)+ (mrar40/mrc)*exp(.446*4*in.zt) + (in.pres_o2/in.pres_co2)*exp(.446*12*in.zt) + (in.pres_co/in.pres_co2)*exp(.446*16*in.zt);

%Scaling factors
sp_yield_n  = 2.4; %Jakosky 1994: escaping N particles per sputtering reaction of N2
sp_yield_c  = .7; %Jakosky 1994: escaping C particles per sputtering reaction of CO2
sp_yield_ar = 1.4; %Jakosky 1994: escaping Ar particles per sputtering reaction of Ar
sp_diff_n2_co2 = exp(.446*16*in.zt); %diffusion of N2 and CO2
sp_diff_ar36_co2 = exp(.446*8*in.zt); %diffusion of Ar36 and CO2
sp_diff_ar38_co2 = exp(.446*6*in.zt); %diffusion of Ar38 and CO2
sp_diff_ar40_co2 = exp(.446*4*in.zt); %diffusion of Ar40 and CO2

%Flux calculations
sput_c = -in.f_sputtering*co2_conv*sput_co2_pps*(1/dilute); %[mbar/myr]

sput_n = -in.f_sputtering*(sp_yield_n/sp_yield_c)*(pn*44/(pc*28))*sp_diff_n2_co2*0.5*n2_conv*sput_co2_pps*(1/dilute); %[mbar/myr] 
%^^^ note the 0.5 is there because "n2_conv" uses M = 28 amu, but the
%nitrogen yield uses particles of N, which have M = 14 amu

sput_ar36 = -in.f_sputtering*(sp_yield_ar/sp_yield_c)*(par36*44/(pc*36))*sp_diff_ar36_co2*ar36_conv*sput_co2_pps*(1/dilute); %[mbar/myr] 
sput_ar38 = -in.f_sputtering*(sp_yield_ar/sp_yield_c)*(par38*44/(pc*38))*sp_diff_ar38_co2*ar38_conv*sput_co2_pps*(1/dilute); %[mbar/myr] 
sput_ar40 = -in.f_sputtering*(sp_yield_ar/sp_yield_c)*(par40*44/(pc*40))*sp_diff_ar40_co2*ar40_conv*sput_co2_pps*(1/dilute); %[mbar/myr] 

%% Outgassing [CNAR]

%Multiply outgassing rate by concentration
outgas_c  = in.og_profile(1,idx)*in.mag_co2;
outgas_n  = in.og_profile(1,idx)*in.mag_n2;

outgas_ar36 = in.og_profile(1,idx)*in.mag_ar36;
outgas_ar38 = in.og_profile(1,idx)*in.mag_ar38;

idx_ar40  = floor((round(sv.t,2) - 700) + 1); 
outgas_ar40 = in.og_profile(1,idx)*(in.ar40m_profile(idx_ar40,1)/in.mantle_mass);

%% Photochemical Loss [CN]

%Euv flux is defined in the sputtering section, for N2
lyman_flux = power(sv.t/4500,-0.86 + in.lym_uncertainty); %Lymam continuum flux, for CO2

%Carbon
photo_c = -in.f_photoC*(lyman_flux^in.lym_pl)*co2_conv*in.photo_flux_c*mrc/in.pres_co2;

%Nitrogen
photoPD_n = -in.f_photoN*(euv_flux^in.euv_pl)*0.5*n2_conv*in.photo_fluxPD_n*mrn/in.pres_n2;
photoDR_n = -in.f_photoN*(euv_flux^in.euv_pl)*0.5*n2_conv*in.photo_fluxDR_n*mrn/in.pres_n2;
photoOP_n = -in.f_photoN*(euv_flux^in.euv_pl)*0.5*n2_conv*in.photo_fluxOP_n*mrn/in.pres_n2;
%^^^ note the 0.5 is there because "n2_conv" uses M = 28 amu

%% Ion Loss [CN]

ionEUV_flux = power(sv.t/4500,-3.51 + 0);

ion_c = -in.f_ion*in.cphi0*ionEUV_flux*co2_conv*mrc/in.pres_co2;
ion_n = -in.f_ion*in.cphi0*ionEUV_flux*0.5*n2_conv*in.ion_ratio*mrn/in.pres_n2;
%^^^ note the 0.5 is there because "n2_conv" uses M = 28 amu

%% Carbonate and Nitrate Deposition [CN]

dep_c = -in.carbonate_profile(1,idx);
dep_n = -in.nitrate_profile(1,idx);

%% Radioactive Decay and Crustal Erosion [Ar K]

% Crustal erosion
erode_ar40 = in.crustal_erosion_rate;

%% Interplanetary Dust Particles [AR]

idp_ar36 = in.idp_rate_36;
idp_ar38 = in.idp_rate_38;

%% Format the function output

c_rates = [sput_c; outgas_c; photo_c; ion_c; dep_c];
n_rates = [sput_n; outgas_n; photoPD_n; photoDR_n; photoOP_n; ion_n; dep_n];
ar36_rates = [sput_ar36; outgas_ar36; idp_ar36];
ar38_rates = [sput_ar38; outgas_ar38; idp_ar38];
ar40_rates = [sput_ar40; outgas_ar40; erode_ar40];

end