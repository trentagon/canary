function [ar40m_new, ar40c_new, km_new, kc_new] = advanceMantle(sv,in,rates)

%Forward Euler advance:
%advance the crust and mantle reservoirs for ar40 and K
%by default all rates for mantle and crust processes are positive

ar40m_new = sv.ar40m + (sum([-rates.ar40(3),rates.ar40(7)])*sv.k);
ar40c_new = sv.ar40c + (sum([-rates.ar40(5),rates.ar40(6)])*sv.k);

km_new = sv.km + (sum([-rates.k(1),-rates.k(3)])*sv.k);
kc_new = sv.kc + (sum([rates.k(1),-rates.k(2)])*sv.k);

end