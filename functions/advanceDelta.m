function [dc_new,dn_new] = advanceDelta(sv,in,rates)

%----Getting the order of rates right is absolutely crucial here!!!!
%Or else, fractionation factors and other delta treatment won't match the
%the rate it is supposed to and the calculation will be wrong

%If the atmosphere is negative this timestep, reset delta to 1
%Also check if a species is held static

%Advance carbon
if in.static.c
    dc_new = sv.dc;
elseif sv.negative_flag.c
    dc_new = 1;
else
    pc = sv.pc;
    dc_new = sv.dc;
    
    diff_c = exp(-.446*in.zt); %Diffusive seperation of C isotopes
    %Negative because 13C is heavier that 12C
    
    %Get the order within the rate vector correct
    sput_idx = 1;
    outgas_idx = 2;
    photo_idx = 3;
    ion_idx = 4;
    dep_idx = 5;
    
    del_pc = abs(rates.c)*sv.k; %change in pressure removed by each process
    
    %Sputtering: rayleigh fractionation
    dc_new = dc_new*power(1 - (del_pc(sput_idx)/pc), diff_c - 1);
    
    %Photochemical: rayleigh fractionation
    dc_new = dc_new*power(1 - (del_pc(photo_idx)/pc), (in.alph_photo_c * diff_c) - 1);
    
    %Ion: rayleigh fractionation
    dc_new = dc_new*power(1 - (del_pc(ion_idx)/pc), diff_c - 1);
    
    %Carbonate: rayleigh fractionation
    dc_new = dc_new*power(1 - (del_pc(dep_idx)/pc),in.alph_carbonate - 1);
    
    %Outgassing: mixing of isotopes
    dc_new = ((del_pc(outgas_idx)*in.dc_mantle)+(pc*dc_new)) / (del_pc(outgas_idx)+pc);
end

%Advance nitrogen
if in.static.n
    dn_new = sv.dn;
elseif sv.negative_flag.n
    dn_new = 1;
else
    pn = sv.pn;
    dn_new = sv.dn;
    
    diff_n = exp(-.446*in.zt); %Diffusive seperation of N isotopes
    %Negative because 15N is heavier that 14N
    
    %Get the order within the rate vector correct
    sput_idx = 1;
    outgas_idx = 2;
    photoPD_idx = 3; %Photodissociation
    photoDR_idx = 4; %Dissociative Recombination
    photoOP_idx = 5; %Other Processes
    ion_idx = 6;
    dep_idx = 7;
    
    del_pn = abs(rates.n)*sv.k; %delta pressure removed by each process

    %Sputtering: rayleigh fractionation
    dn_new = dn_new*power(1 - (del_pn(sput_idx)/pn), diff_n - 1);
    
    %Photodissociation: rayleigh fractionation
    dn_new = dn_new*power(1 - (del_pn(photoPD_idx)/pn), (in.alph_photoPD_n * diff_n) - 1);
    
    %Dissociative Recombination: rayleigh fractionation
    dn_new = dn_new*power(1 - (del_pn(photoDR_idx)/pn), (in.alph_photoDR_n * diff_n) - 1);
    
    %Other Processes: rayleigh fractionation
    dn_new = dn_new*power(1 - (del_pn(photoOP_idx)/pn), (in.alph_photoOP_n * diff_n) - 1);
    
    %Ion: rayleigh fractionation
    dn_new = dn_new*power(1 - (del_pn(ion_idx)/pn), diff_n - 1);
    
    %Carbonate: rayleigh fractionation
    dn_new = dn_new*power(1 - (del_pn(dep_idx)/pn),in.alph_nitrate - 1);
    
    %Outgassing: mixing of isotopes
    dn_new = ((del_pn(outgas_idx)*in.dn_mantle)+(pn*dn_new)) / (del_pn(outgas_idx)+pn);    
end

%Take real part
dc_new = real(dc_new);
dn_new = real(dn_new);

end