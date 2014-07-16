function [h] = shorthairH2Oadjust(time, h2o)

%The irga at shorthair has a bad calibration for several periods
%We are going to adjust the irga h2o calibration before we calculate the fluxes

%There are 2 periods
%1) time < 461
%the quadratic function that fits x=irga h2o, y=hmp h2o:
%hmp h2o = 0.03508 * (irga h2o)^2 + 0.202 * (irga h2o) + 2.235

%2) time > 461
%the quadratic function that fits x=irga h2o, y=hmp h2o:
%hmp h2o = 0.0056104* (irga h2o)^2 + 0.08938 * (irga h2o) + 0.53936


%input
%time = current interval time
%h2o = irga h2o mmol/mol moist air

%output
%h = h2o mmol/mol moist air adjusted to hmp h2o


%shorthair h2o correction lookup

%first half correction
if time > 339 && time < 354;
    h=(0.003508 * h2o .* h2o) + (0.202 * h2o) + 2.235;
elseif time > 372.45 && time < 374.42;
    h=(0.003508 * h2o .* h2o) + (0.202 * h2o) + 2.235;
elseif time > 379 && time < 391.42;
    h=(0.003508 * h2o .* h2o) + (0.202 * h2o) + 2.235;
elseif time > 399.8 && time < 407;
    h=(0.003508 * h2o .* h2o) + (0.202 * h2o) + 2.235;
elseif time > 435 && time < 461;
    h=(0.003508 * h2o .* h2o) + (0.202 * h2o) + 2.235;
%second half of the time series
elseif time > 679.7 && time < 685.82;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 695 && time < 696.38;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 700 && time < 702.38;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 710.4 && time < 711.4;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 716.5 && time < 717.4;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 763.4 && time < 764.38;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 783.4 && time < 784.39;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 790.5 && time < 791.5;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 793.6 && time < 794.38;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 797.6 && time < 798.38;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 881.4 && time < 882.4;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 971.36 && time < 973.4;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 1018.8 && time < 1019.38;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 1044.5 && time < 1045.38;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 1048.4 && time < 1049.4;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 1099.4 && time < 1100.4;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 1105.86 && time < 1106.4;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
elseif time > 1120.45 && time < 1121.4;
    h= (0.0056104 * h2o .* h2o) + (0.08938 * h2o) + 0.53936;
else
    h=h2o;
end