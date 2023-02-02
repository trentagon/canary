
t_start = 100;
t_end = 4500;

og_totals = 10^6 * [74,168,148,150,81,20]; % ext:int = 8.5:1, from Greeley and Schneid 2001
chron_ivan = t_end - 1000*[3.86,3.74,3.65,3.46,1.45,0.387];
chron_hart = t_end - 1000*[3.85,3.57,3.40,3.00,0.880,0.235];

chron = chron_hart;

og_smooth = 0.2;

og_profile = smooth_profile(og_totals,chron,t_start,t_end,og_smooth); %km^3/Myr
og_profile = grottify(og_profile,t_start);

semilogy(og_profile(2,:),og_profile(1,:)*10^-6);