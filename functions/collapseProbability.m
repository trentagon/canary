function fc = collapseProbability(sv, in)
%{
This function returns the probability that the atmosphere is collapsed when
taking into account obliquity variations from Laskar and the critical
pressure as a function of obliquity from Forget
%}

maxp = max(in.obliquity_data.c_i(:,1)); %Maximum pressure collapse can occur at
ptot = (sv.pc + sv.pn + sv.par36 + sv.par38 + sv.par40)/1000; %Model total pressure

if ptot > maxp
    %If the current pressure is greater than the maximum collapse pressure
    fc = 0; % 0% chance the atmosphere collapses, regardless of obliquity
else
    
    %Pick guassian from Laskar that is closest to the model time
    t_a = in.t_end - sv.t;
    t_difs = abs(in.obliquity_data.profile.times - t_a);
    t_idx = find(t_difs == min(t_difs));
    if length(t_idx) > 1
        t_idx = t_idx(end);
    end
    
    %Load the gaussian
    y = in.obliquity_data.profile.pdfs(t_idx,:);
    x = in.obliquity_data.profile.x;
    
    %Find closest pressure value from Forget
    p_difs = abs(in.obliquity_data.c_i(:,1) - ptot);
    p_idx = find(p_difs == min(p_difs));
    if length(p_idx) > 1
        p_idx = p_idx(end);
    end
    
    %Get the critical obliquity
    c_i = in.obliquity_data.c_i(p_idx,2); 
    % i.e. below this obliquity, the atmosphere is collapsed
    
    %Find the critical obliquity in the laskar gaussians
    int_cap = find(abs(x-c_i) == min(abs(x-c_i))); 
    if length(int_cap) > 1
        int_cap = int_cap(end);
    end
    
    %Integrate the gaussian pdf up to the critical obliquity
    fc = trapz(x(1:int_cap),y(1:int_cap)); 
    %This is the probability that obliquity was less than critical
    % i.e. probability that the atmosphere is collapsed
    
end

end