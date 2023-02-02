function canary_plot(out)

tmyr = out.times.t;
tstep_myr = out.times.k;
tga = out.times_Ga.t;
tstep_ga = out.times_Ga.k;

pco2 = out.pressures.c;
pn2 = out.pressures.n;
par36 = out.pressures.ar36;
par38 = out.pressures.ar38;
par40 = out.pressures.ar40;
partot = out.par_total;
ptot = pco2 + pn2 + out.par_total;

dco2 = out.deltas_pm.c;
dn2 = out.deltas_pm.n;
dar38 = out.deltas_pm.ar38;
dar40 = out.deltas_pm.ar40;

flux_co2 = abs(out.rates.c);
flux_n2 = abs(out.rates.n);

ar_sput = out.rates.ar36(:,1)+out.rates.ar38(:,1)+out.rates.ar40(:,1);
ar_og_36_38 = out.rates.ar36(:,2)+out.rates.ar38(:,2);
ar_og_40 = out.rates.ar40(:,2);
ar_ce = out.rates.ar40(:,3);
ar_idp = out.rates.ar36(:,3) + out.rates.ar38(:,3);
flux_ar = abs([ar_sput,ar_og_36_38,ar_og_40,ar_ce, ar_idp]);

fsum_co2 = round(sum(tstep_myr.*flux_co2,1),2);
fsum_n2 = round(sum(tstep_myr.*flux_n2,1),2);
fsum_ar = sum(tstep_myr.*flux_ar,1);

co2_labels = ["Sputtering","Outgassing","Photochemical (all)","Ion","Carbonate Deposition"];
n2_labels = ["Sputtering","Outgassing","Photochemical (Photodissociation)","Photochemical (Dissociative Recombination)","Photochemical (Other Processes)","Ion","Carbonate Deposition"];
ar_labels = ["Sputtering","Outgassing","Crustal Erosion"];

subplot = @(m,n,p) subtightplot (m, n, p, [0.075 0.065], [0.1 0.01], [0.1 0.01]);
fg = figure('Renderer', 'painters', 'Position', [10 10 1800 1000]);

blue1 = [158, 194, 255]/255;
blue2 = [68, 117, 201]/255;
blue3 = [7, 59, 148]/255;
green1 = [144, 230, 133]/255;
green2 = [54, 138, 44]/255;
green3 = [14, 71, 7]/255;
black1 = [0,0,0];
red1 = [161, 63, 63]/255;
orange1 = [214, 159, 64]/255;

ca = [86,180,233]/255;
cb = [0,158,115]/255;
cc = [240,228,66]/255;
cd = [204,121,167]/255;
ce = [213,94,0]/255;

co2_color = black1;
n2_color = blue2;
ar_color = red1;
ar_color2 = green2;

color_sput = ce;
color_photo = cb;
color_ion = cc;
color_dep = cd;
color_og = black1;
color_crust = ca;

font = 26;
ticks_font_size = 22;
legend_font_size = 18;
line_width = 3;

subplot(2,3,1);
semilogy(tga,pco2,'LineWidth', line_width, 'Color', co2_color, 'DisplayName','pCO_{2}');
hold on
semilogy(tga,pn2,'LineWidth', line_width, 'Color', n2_color, 'DisplayName','pN_{2}');
semilogy(tga,partot,'LineWidth', line_width, 'Color', ar_color, 'DisplayName','pAr');
%semilogy(tga,ptot,'LineWidth', line_width, 'Color', black, 'DisplayName','pTOT');
xlabel('Time before present (Ga)','FontSize',font);
ylabel('Partial pressure (mbar)','FontSize',font);
set(gca,'Xdir','reverse');
set(gca,'FontSize',ticks_font_size,'XMinorTick','on','YMinorTick','on');
set(gca,'xtick',[0,.5:0.5:3.5,3.8],'Xlim',[0,3.8],'Ylim',[10^(-4),10^4]);
box on
ew = 3;
errorbar(-0.04, 34, 20, 'Color', co2_color, 'LineWidth', ew, 'DisplayName','','clipping','off'); %C
errorbar(-0.04, 0.25, 0.13, 'Color', n2_color, 'LineWidth', ew, 'DisplayName','','clipping','off'); %N
errorbar(-0.08, 0.26, 0.15, 'Color', ar_color, 'LineWidth', ew, 'DisplayName','','clipping','off'); %Ar
legend({'pCO_{2}','pN_{2}','pAr'},'FontSize',legend_font_size);
hold off

subplot(2,3,2);
plot(tga,dco2,'LineWidth', line_width, 'Color', co2_color, 'DisplayName','\deltaC^{13}');
hold on
plot(tga,dn2,'LineWidth', line_width, 'Color', n2_color, 'DisplayName','\deltaN^{15}');
plot(tga,dar38,'LineWidth', line_width, 'Color', ar_color, 'DisplayName','\deltaAr^{38}');
plot(tga,dar40,'LineWidth', line_width, 'Color', ar_color2, 'DisplayName','\deltaAr^{40}');
xlabel('Time before present (Ga)','FontSize',font);
ylabel('Isotopic Composition (per mil)','FontSize',font);
set(gca,'Xdir','reverse');
set(gca,'FontSize',ticks_font_size,'XMinorTick','on','YMinorTick','on');
set(gca,'xtick',[0,.5:0.5:3.5,3.8],'Xlim',[0,3.8],'Ylim',[-1000,6000]);
box on
%Modern values
errorbar(-0.04, 46, 2*3, 'Color', co2_color, 'LineWidth', ew, 'DisplayName','','clipping','off'); %C
errorbar(-0.08, 572, 82*3, 'Color', n2_color, 'LineWidth', ew, 'DisplayName','','clipping','off'); %N
errorbar(-0.04, 310, 30*3, 'Color', ar_color, 'LineWidth', ew, 'DisplayName','','clipping','off'); %Ar38
errorbar(-0.12, 0, 79*3, 'Color', ar_color2, 'LineWidth', ew, 'DisplayName','','clipping','off'); %Ar40
%Ancient values
% errorbar(3.82, 18.475, 18.475, 'Color', ar_color, 'LineWidth', ew, 'DisplayName','','clipping','off'); %Ar38
% errorbar(3.82, -670.53, 52.64, 'Color', ar_color2, 'LineWidth', ew, 'DisplayName','','clipping','off'); %Ar40
legend({'\delta^{13}C','\delta^{15}N','\delta^{38}Ar','\delta^{40}Ar'},'FontSize',legend_font_size);
hold off

subplot(2,3,4);
semilogy(tga,flux_co2(:,1),'LineWidth', line_width, 'Color', color_sput, 'DisplayName',['Sputtering: ',num2str(fsum_co2(1)),' mbar']);
hold on
semilogy(tga,flux_co2(:,3),'LineWidth', line_width, 'Color', color_photo, 'DisplayName',['Photochemical: ',num2str(fsum_co2(3)),' mbar']');
semilogy(tga,flux_co2(:,4),'LineWidth', line_width, 'Color', color_ion, 'DisplayName',['Ion: ',num2str(fsum_co2(4)),' mbar']');
semilogy(tga,flux_co2(:,5),'LineWidth', line_width, 'Color', color_dep, 'DisplayName',['Carbonate Deposition: ',num2str(fsum_co2(5)),' mbar']');
semilogy(tga,flux_co2(:,2),'LineWidth', line_width, 'Color', color_og, 'DisplayName',['Outgassing: ',num2str(fsum_co2(2)),' mbar']');
xlabel('Time before present (Ga)','FontSize',font);
ylabel('pCO_{2} Flux (mbar Myr^{-1})','FontSize',font);
set(gca,'Xdir','reverse');
set(gca,'FontSize',ticks_font_size,'XMinorTick','on','YMinorTick','on');
set(gca,'xtick',[0,.5:0.5:3.5,3.8],'ytick',[10^-10,10^-8,10^-6,10^-4,10^-2,10^0,10^2],'Xlim',[0,3.8],'Ylim',[10^-10,10^2]);
box on
legend({},'FontSize',legend_font_size);
hold off

subplot(2,3,5);
semilogy(tga,flux_n2(:,1),'LineWidth', line_width, 'Color', color_sput, 'DisplayName',['Sputtering: ',num2str(fsum_n2(1)),' mbar']);
hold on
semilogy(tga,flux_n2(:,3)+flux_n2(:,4)+flux_n2(:,5),'LineWidth', line_width, 'Color', color_photo, 'DisplayName',['Photochemical: ',num2str(round(fsum_n2(3)+fsum_n2(4)+fsum_n2(5),2)),' mbar']);
semilogy(tga,flux_n2(:,6),'LineWidth', line_width, 'Color', color_ion, 'DisplayName',['Ion: ',num2str(fsum_n2(6)),' mbar']);
semilogy(tga,flux_n2(:,7),'LineWidth', line_width, 'Color', color_dep, 'DisplayName',['Nitrate Deposition: ',num2str(fsum_n2(7)),' mbar']);
semilogy(tga,flux_n2(:,2),'LineWidth', line_width, 'Color', color_og, 'DisplayName',['Outgassing: ',num2str(fsum_n2(2)),' mbar']);
xlabel('Time before present (Ga)','FontSize',font);
ylabel('pN_{2} Flux (mbar Myr^{-1})','FontSize',font);
set(gca,'Xdir','reverse');
set(gca,'FontSize',ticks_font_size,'XMinorTick','on','YMinorTick','on');
set(gca,'xtick',[0,.5:0.5:3.5,3.8],'ytick',[10^-10,10^-8,10^-6,10^-4,10^-2,10^0,10^2],'Xlim',[0,3.8],'Ylim',[10^-10,10^2]);
box on
legend({},'FontSize',legend_font_size);
hold off

subplot(2,3,6);
semilogy(tga,flux_ar(:,1),'LineWidth', line_width, 'Color', color_sput, 'DisplayName',['Sputtering: ',num2str(fsum_ar(1),'%.2g'),' mbar']);
hold on
semilogy(tga,flux_ar(:,4),'LineWidth', line_width, 'Color', color_crust, 'DisplayName',['Crustal Erosion: ',num2str(fsum_ar(4),'%.2g'),' mbar']);
semilogy(tga,flux_ar(:,5),'LineWidth', line_width, 'Color', color_ion, 'DisplayName',['IDP Accretion: ',num2str(fsum_ar(5),'%.2g'),' mbar']);
semilogy(tga,flux_ar(:,2),'LineWidth', line_width, 'Color', color_og, 'DisplayName',['^{36}Ar + ^{38}Ar Outgassing: ',num2str(fsum_ar(2),'%.2g'),' mbar']);
semilogy(tga,flux_ar(:,3),'LineWidth', line_width, 'Color', color_dep, 'DisplayName',['^{40}Ar Outgassing: ',num2str(fsum_ar(3),'%.2g'),' mbar']);
xlabel('Time before present (Ga)','FontSize',font);
ylabel('pAr Flux (mbar Myr^{-1})','FontSize',font);
set(gca,'Xdir','reverse');
set(gca,'FontSize',ticks_font_size,'XMinorTick','on','YMinorTick','on');
set(gca,'xtick',[0,.5:0.5:3.5,3.8],'ytick',[10^-10,10^-8,10^-6,10^-4,10^-2,10^0,10^2],'Xlim',[0,3.8],'Ylim',[10^-10,10^2]);
box on
legend({},'FontSize',legend_font_size);
hold off


repcase_parms = [out.inputs.f_outgassing, out.inputs.f_sputtering, out.inputs.f_photoC, out.inputs.f_photoN, out.inputs.f_ce, out.inputs.zt, out.inputs.mag_co2, out.inputs.euv_pl,out.inputs.lym_pl];
repcase_parms(7) = repcase_parms(7)*10^6;
repcase_params = round(repcase_parms,2);

par_labels = {'f_{og}','f_{sp}','f_{pC}','f_{pN}','f_{ce}','z/T','X^{m}_{CO2}','euv_{pl}','lym_{pl}'};

dim = [.732 .68 .3 .3];
param_cap = {['f_{og} = ', num2str(repcase_params(1))],['f_{sp} = ', num2str(repcase_params(2))], ['f_{pC} = ', num2str(repcase_params(3))], ['f_{pN} = ', num2str(repcase_params(4))],['f_{ce} = ', num2str(round(out.inputs.f_ce,3))], ['\Deltaz/T = ', num2str(repcase_params(6)), ' km K^{-1}'], ['X^{C}_{mag} = ', num2str(repcase_params(7)), ' ppm'], ['a_{euv} = ', num2str(repcase_params(8))], ['a_{lym} = ', num2str(repcase_params(9))]};
annotation('textbox',dim,'String',param_cap,'FitBoxToText','on', 'FontSize',font);

AddLetters2Plots(fg)


end