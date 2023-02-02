function [ar40m_profile, crustal_erosion_rate] = calculate_40Ar_profiles(in)

%Load in full outgassing profile from 4.4 to 0 Ga depending on chronology,
%and outgassing multiplier

og_profile_ar_in_function(1,:) = in.og_profile_ar(1,:)*in.f_outgassing*in.crust_density*in.og_convert; % [kg/Myr] 1e9 converts density from m^-3 to km^-3
og_profile_ar_in_function(2,:) = in.og_profile_ar(2,:);

%Evolve system forward through full evolution to get profiles

t_pre = 100; % start at 4.4 Ga = 100 Myr in the model timekeeping
k = 1; % timestep = 1 Myr for precalculation

K_44ga = 0.4e-6; % K40 concentration in mantle at 4.4 Ga

Km0 = in.mantle_mass*K_44ga; % initial K in mantle
Kc0 = 0; % initial K in crust

Arm0 = 0; % initial Ar40 in mantle
Arc0 = 0; % initial Ar40 in crust

[sv_pre.t, sv_pre.Km, sv_pre.Kc, sv_pre.Arm, sv_pre.Arc] = deal(t_pre, Km0, Kc0, Arm0, Arc0);

ar40m_rec = [];
t_rec = [];

%Loop to integrate Ar40 parameters from 4.4 to 3.8 Ga
running = 1;
while running
    
    % Calculate rates
    
    idx  = floor(100*(round(sv_pre.t,2) - t_pre) + 1); % get discrete idx closest to model time
    emp_rate = og_profile_ar_in_function(1,idx)/(in.mars_g * 0.01 / (4 * pi * (in.mars_r)^2)); % volcanic emplacement rate [kg/Myr]
        
    Ar40_og_kg = (sv_pre.Arm/in.mantle_mass) * emp_rate * k; %  kg Ar / Myr
    K_og_kg = (5*sv_pre.Km/in.mantle_mass) * emp_rate* k; % kg K/kg * kg/Myr = kg K / Myr
    
    new_Km = sv_pre.Km + k*(-sv_pre.Km*in.lambda_k40 - K_og_kg); % Update K in mantle
    new_Kc = sv_pre.Kc + k*(-sv_pre.Kc*in.lambda_k40 + K_og_kg); % Update K in crust
        
    new_Arm = sv_pre.Arm + k*(in.lambda_k40*in.beta_ar40*sv_pre.Km - Ar40_og_kg); % Update Ar40 in mantle
    new_Arc = sv_pre.Arc + k*(in.lambda_k40*in.beta_ar40*sv_pre.Kc); % Update Ar40 in crust
    
    new_t = sv_pre.t + k; % Update time
    
    [sv_pre.t, sv_pre.Km, sv_pre.Kc, sv_pre.Arm, sv_pre.Arc] = deal(new_t, new_Km, new_Kc, new_Arm, new_Arc);
    
    ar40m_rec(end+1,1) = sv_pre.Arm;
    t_rec(end+1,1) = sv_pre.t;
    
    % Check for end of run
    if sv_pre.t >= 4500
        break
    end
    
end

idx = (in.t_start - 100)/k;

ar40m_profile = [ar40m_rec(idx:end),t_rec(idx:end)]; %kg
crustal_erosion_rate = (in.f_ce*sv_pre.Arc/(4400)) * in.mars_g * 0.01 / (4 * pi * (in.mars_r)^2); %mbar
% crustal_erosion_rate = (in.f_ce*sv_pre.Arc/(in.t_end-in.t_start)); %kg

end