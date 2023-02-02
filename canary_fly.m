%{

    _____          _   _          _______     __
   / ____|   /\   | \ | |   /\   |  __ \ \   / /   \\
  | |       /  \  |  \| |  /  \  | |__) \ \_/ /    (o>
  | |      / /\ \ |   ` | / /\ \ |  _  / \   /  \\_//)
  | |____ / ____ \| |\  |/ ____ \| | \ \  | |    \_/_)
   \_____/_/    \_\_| \_/_/    \_\_|  \_\ |_|     _|_

A Mars Atmospheric Evolution Model that tracks [CA]rbon, [N]itrogen, and [AR]gon                                         
Developed by Trent Thomas and Renyu Hu 

%}
%% Model Startup

%Clear everything
close all
clear all

%Set Path
restoredefaultpath
canary_path = '/Users/tthomas/Documents/research/mars-canary/CANARY/'; %Path to the CANARY directory
addpath(genpath(canary_path)); %Add all subdirectories to path

%Start clock
tStart = tic;

%% Load in input values:
run('inputFile_fig5.m'); %Run input file
clearvars -except in tStart; %Only keep the input struct and time

%% Initialize Data Structures

%----Core structs----
%"in": inputs, contains science inputs and a few technical specifications
%"sv": state vector, contains primary variables that update during each loop
%"rates": rates, contains sources and sinks for each species
%"out": output, contains saved values for outputting

rates = struct();

%Initialize output
out.title = "This is the output data from a run of canary_main.m";
out.runDate = datestr(now, 'mm/dd/yy-HH:MM:SS');

%Create storage
%Atmosphere
out.pressures.c = []; 
out.pressures.n = []; 
out.pressures.ar36 = [];
out.pressures.ar38 = [];
out.pressures.ar40 = [];
out.deltas.c = []; 
out.deltas.n = [];
out.deltas.ar38 = [];
out.deltas.ar40 = [];
%Rates
out.rates.c  = []; 
out.rates.n  = []; 
out.rates.ar36 = [];
out.rates.ar38 = []; 
out.rates.ar40 = [];
%Times
out.times.t = []; 
out.times.k = [];
%Negative atmospheres and collapse
out.negative_flag.c = [];
out.negative_flag.n = [];
out.negative_flag.ar36 = [];
out.negative_flag.ar38 = [];
out.negative_flag.ar40 = [];
out.collapse = [];

%% Backward integrate for initial pressures

%Initialize state vector
[sv.pc, sv.pn, sv.par36, sv.par38, sv.par40] = deal(34, 0.25, in.min_pressure_ar, in.min_pressure_ar, 0.2598); %Initial pressures
sv.t = in.t_end; %time
sv.fc = 0; %This will change if atmospheric collapse is treated normally

%Ignore 36Ar and 38Ar during backward integration
rates.ar36 = 0;
rates.ar38 = 0;

running_backward = 1;
while running_backward
        
    %%%%%%%%%%%%%%%%%%%%%% RATE CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%
    %Calculate rate of impact of the sources and sinks
    %The switch goes through the different methods of handling collapse
    switch in.collapse
        case "obliquity"
            [rates.c, rates.n, ~, ~, rates.ar40] = calculateRates(sv.pc, sv, in); %Uncollapsed atmosphere [pCO2 = sv.pc]
            sv.fc = collapseProbability(sv,in); %Calculate the probability that the atmosphere is collapsed
            
            %Is there a chance the atmosphere is collapsed?
            if sv.fc > 0
                [col_c, col_n, ~, ~, col_ar40] = calculateRates(in.pc_col, sv, in); %Collapsed atmosphere [pCO2 = in.pc_col]
                
                %Weighted average of the collapsed and uncollapsed atmospheres
                rates.c = rates.c*(1-sv.fc) + col_c*sv.fc;
                rates.n = rates.n*(1-sv.fc) + col_n*sv.fc;
                rates.ar40 = rates.ar40*(1-sv.fc) + col_ar40*sv.fc;
            end
        case "static"
        
            %Is the total pressure below a set threshold value?
            if sum([sv.pc,sv.pn,sv.par36,sv.par38,sv.par40]) < in.static_collapse_threshold
                [rates.c, rates.n, ~, ~, rates.ar40] = calculateRates(in.pc_col, sv, in); %Collapsed atmosphere [pCO2 = in.pc_col]
            else
                [rates.c, rates.n, ~, ~, rates.ar40] = calculateRates(sv.pc, sv, in); %Uncollapsed atmosphere [pCO2 = sv.pc]
            end
        case "nothing"
            [rates.c, rates.n, ~, ~, rates.ar40] = calculateRates(sv.pc, sv, in); %Uncollapsed atmosphere [pCO2 = sv.pc]
        otherwise
            error("Error: invalid value of in.collapse - must be 'obliquity', 'static', or 'nothing'.");
    end
    
    %%%%%%%%%%%%%%%%%%%%%% TIMESTEP CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%
    %Determine the maximum allowable timestep
    sv.k = calculateTimestep(sv,in,rates);

    %%%%%%%%%%%%%%%%%%%%%% ADVANCE THE STATE VECTOR %%%%%%%%%%%%%%%%%%%%%%
    %Advance the presssure and delta values according to rates and tstep
    
    %Record current pressures
    old_pc = sv.pc;
    old_pn = sv.pn;
    old_par40 = sv.par40;

    %Going backward so rates are negative
    rates.c = -1*rates.c;
    rates.n = -1*rates.n;
    rates.ar40 = -1*rates.ar40;
    
    %Calculate new pressures and deltas based on current state and rates
    [new_pc, new_pn, ~, ~, new_par40, back_flag] = advancePressure(sv,in,rates);
        
    %Update state vector
    [sv.pc, sv.pn, sv.par40] = deal(new_pc, new_pn, new_par40);
    sv.t = sv.t - sv.k;
        
    %%%%%%%%%%%%%%%%%%%%%% CHECK IF RUN IS DONE %%%%%%%%%%%%%%%%%%%%%%
    if sv.t < in.t_start
        running_backward = 0;
    elseif any([back_flag.c;back_flag.n;back_flag.ar40])
        error("Error: negative pressure encountered during backward integration");
    end
    
end

%Calculate 36Ar and 38Ar from mantle and ALH values
old_par36 = old_par40/in.ar40to36_start;
old_par38 = old_par36/in.ar36to38_start;

fprintf("Backward integration results:\n");
fprintf("pCO2 = %.2f | pN2  = %.2f [mbar]\n",old_pc,old_pn);
fprintf("pAr36 = %.2g | pAr38 = %.2g | pAr40 = %.2f [mbar]\n",old_par36,old_par38,old_par40);

%% Run the model loop starting from the backward integration results

%Initialize state vector

dar38_back = (old_par38/old_par36)/(1/5.501);
dar40_back = (old_par40/old_par36)/(1900);

[sv.pc, sv.pn, sv.par36, sv.par38, sv.par40] = deal(old_pc, old_pn, old_par36, old_par38, old_par40); %Initial pressures
[sv.dc, sv.dn, sv.dar38, sv.dar40] = deal(in.dc_0, in.dn_0, dar38_back, dar40_back); %Initial deltas
sv.t = in.t_start; %time
%sv.k is timestep, calculated in first loop
%sv.fc is probability of collapse, calculated in first loop - if necessary
sv.negative_flag.c = 0; sv.negative_flag.n = 0; 
sv.negative_flag.ar36 = 0; sv.negative_flag.ar38 = 0; sv.negative_flag.ar40 = 0;
%negative_flag.X is 1 if the atmosphere goes negative

tLoop = tic; %Record the loop time
running = 1;
while running
    
    %%%%%%%%%%%%%%%%%%%%%% COMET CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%
    %Does a comet reset the atmospheric Argon?
    if in.comet == "yes" && sv.t >= in.impact_time && in.impacted == 1
        sv.par = in.par_comet; %Reset Pressure
        sv.dar = in.dar_comet; %Reset Delta
        in.impacted = 0; %It already impacted, don't do it again
        fprintf("Comet impacted at %.2f Myrs\n",sv.t);
    end
    
    %%%%%%%%%%%%%%%%%%%%%% RATE CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%
    %Calculate rate of impact of the sources and sinks
    %The switch goes through the different methods of handling collapse
    switch in.collapse
        case "obliquity"
            [rates.c, rates.n, rates.ar36, rates.ar38, rates.ar40] = calculateRates(sv.pc, sv, in); %Uncollapsed atmosphere [pCO2 = sv.pc]
            sv.fc = collapseProbability(sv,in); %Calculate the probability that the atmosphere is collapsed
            
            %Is there a chance the atmosphere is collapsed?
            if sv.fc > 0
                [col_c, col_n, col_ar36, col_ar38, col_ar40] = calculateRates(in.pc_col, sv, in); %Collapsed atmosphere [pCO2 = in.pc_col]
                
                %Weighted average of the collapsed and uncollapsed atmospheres
                rates.c = rates.c*(1-sv.fc) + col_c*sv.fc;
                rates.n = rates.n*(1-sv.fc) + col_n*sv.fc;
                rates.ar36 = rates.ar36*(1-sv.fc) + col_ar36*sv.fc;
                rates.ar38 = rates.ar38*(1-sv.fc) + col_ar38*sv.fc;
                rates.ar40 = rates.ar40*(1-sv.fc) + col_ar40*sv.fc;
            end
        case "static"
        
            %Is the total pressure below a set threshold value?
            if sum([sv.pc,sv.pn,sv.par36,sv.par38,sv.par40]) < in.static_collapse_threshold
                [rates.c, rates.n, rates.ar36, rates.ar38, rates.ar40] = calculateRates(in.pc_col, sv, in); %Collapsed atmosphere [pCO2 = in.pc_col]
            else
                [rates.c, rates.n, rates.ar36, rates.ar38, rates.ar40] = calculateRates(sv.pc, sv, in); %Uncollapsed atmosphere [pCO2 = sv.pc]
            end
        case "nothing"
            [rates.c, rates.n, rates.ar36, rates.ar38, rates.ar40] = calculateRates(sv.pc, sv, in); %Uncollapsed atmosphere [pCO2 = sv.pc]
        otherwise
            error("Error: invalid value of in.collapse - must be 'obliquity', 'static', or 'nothing'.");
    end
        
    %%%%%%%%%%%%%%%%%%%%%% TIMESTEP CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%
    %Determine the maximum allowable timestep
    sv.k = calculateTimestep(sv,in,rates);
    
    %%%%%%%%%%%%%%%%%%%%%% WRITE TO OUTPUT STRUCT %%%%%%%%%%%%%%%%%%%%%%
    %This function is expensive because of rate recording but preallocation
    %doesn't fix it so I don't know how to make this better
    out = writeOutput(sv,rates,out);

    %%%%%%%%%%%%%%%%%%%%%% ADVANCE THE STATE VECTOR %%%%%%%%%%%%%%%%%%%%%%
    %Advance the presssure and delta values according to rates and tstep
    
    %Calculate new pressures and deltas based on current state and rates
    [new_pc, new_pn, new_par36, new_par38, new_par40, sv.negative_flag] = advancePressure(sv,in,rates);
    [new_dc, new_dn] = advanceDelta(sv,in,rates);
    
    %Manually calculate deltas for Ar. Relative to solar and modern Mars.
    new_dar38 = (new_par38/new_par36)/(1/5.501);
    new_dar40 = (new_par40/new_par36)/(1900);
    
    %Update state vector
    [sv.pc, sv.pn, sv.par36, sv.par38, sv.par40] = deal(new_pc, new_pn, new_par36, new_par38, new_par40);
    [sv.dc, sv.dn, sv.dar38, sv.dar40] = deal(new_dc, new_dn, new_dar38, new_dar40);
    sv.t = sv.t + sv.k;
        
    %%%%%%%%%%%%%%%%%%%%%% CHECK IF RUN IS DONE %%%%%%%%%%%%%%%%%%%%%%
    running = checkHalt(sv,in);
    
end

%% Process model output
out.pressures.units = "mbar";
out.deltas.units = "unitless [model isotope ratio/standard isotope ratio]";
out.rates.units = "mbar per Myr";
out.rates.c_labels = ["Sputtering";"Outgassing";"Photochemical (all)";"Ion";"Carbonate Deposition"]';
out.rates.n_labels = ["Sputtering";"Outgassing";"Photochemical (Photodissociation)";"Photochemical (Dissociative Recombination)";"Photochemical (Other Processes)";"Ion";"Carbonate Deposition"]';
out.rates.ar_labels = ["Sputtering";"Outgassing";"Crustal Erosion"]';
out.rates.ar36_labels = ["Sputtering 36";"Outgassing 36"];
out.rates.ar38_labels = ["Sputtering 38";"Outgassing 38"];
out.rates.ar40_labels = ["Sputtering 40"; "Outgassing 40 [mbar]"; "Crust Erosion 40 [mbar]"];
out.times.units = "Myr [million years]";
out.times.what_is_k = "k is the timestep. it has nothing to do with potassium. same units as time, friend.";
out.inputs = in;
runtime_total = toc(tStart); %total runtime
runtime_loop  = toc(tLoop); %runtime inside the main loop
out.runtime_total = [num2str(runtime_total), ' seconds'];
out.runtime_loop = [num2str(runtime_loop), ' seconds'];

%% Make secondary calculations
out.title2 = "\/--- secondary calculations/conversions below ---\/";

%Delta values in per mil
out.deltas_pm.c = (out.deltas.c - 1)*1000;
out.deltas_pm.n = (out.deltas.n - 1)*1000;
out.deltas_pm.ar38 = (out.deltas.ar38 - 1)*1000;
out.deltas_pm.ar40 = (out.deltas.ar40 - 1)*1000;
out.deltas_pm.units = "per mil";

%Did the run ever go negative for any species?
out.any_negative = [any(out.negative_flag.c);any(out.negative_flag.n);any(out.negative_flag.ar36);any(out.negative_flag.ar38);any(out.negative_flag.ar40)];
if any(out.any_negative)
    fprintf("Warning: atmosphere went negative during model run");
end

%Convert Myr to Ga
out.times_Ga.t = (in.t_end - out.times.t)/1000;
out.times_Ga.k = out.times.k/1000;
out.times_Ga.units = "Ga [billion years ago]";

%Sum argon
out.par_total = out.pressures.ar36 + out.pressures.ar38 + out.pressures.ar40;

%% Summary, Plots, and Analysis

%Runtime
fprintf("\nRuntime = %.2f seconds\n", runtime_total);

%Report the model results
fprintf("------------------------------------------------------------------\n");
fprintf("~~** Model Results **~~\n")
fprintf("CARBON:     p_start = %.4f | p_finish = %.4f | delta = %.4f\n", old_pc, out.pressures.c(end), out.deltas_pm.c(end));
fprintf("NITROGEN:   p_start = %.4f | p_finish = %.4f | delta = %.4f\n", old_pn, out.pressures.n(end), out.deltas_pm.n(end));
fprintf("ARGON 36:   p_start = %.4g | p_finish = %.4g\n", old_par36, out.pressures.ar36(end));
fprintf("ARGON 38:   p_start = %.4g | p_finish = %.4g\n", old_par38, out.pressures.ar38(end));
fprintf("ARGON 40:   p_start = %.4f | p_finish = %.4f\n", old_par40, out.pressures.ar40(end));
fprintf("ARGON TOT:  p_start = %.4f | p_finish = %.4f | delta38 = %.4f | delta40 = %.4f\n", out.par_total(1), out.par_total(end), out.deltas_pm.ar38(end), out.deltas_pm.ar40(end));

%Print target values
fprintf("------------------------------------------------------------------\n");
fprintf(":+: Present Mars [1 sigma] :+:\n");
fprintf("pCO2 = %.2f +/- %.3f | d13C  = %.2f  +/- %.f\n",34,6.67,46,4);
fprintf("pN2  = %.2f  +/- %.3f | d15N  = %.2f +/- %.f\n",0.25,0.043,572,82);
fprintf("pAr  = %.2f  +/- %.3f | d38Ar = %.2f +/- %.f | d40Ar = %.2f +/- %.f\n\n",0.26,0.05,310,30,0,157.9);

%Summary plots
canary_plot(out);
