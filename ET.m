%Inputs are eddy covariance data at G:\Aaron\Final\ECdata\
%------------------------------------------------------------------------%
%Outputs
% ET gives the cumulative ET in mm using regress and robust functions
%------------------------------------------------------------------------%
clear all

disp('Work on this site:')
%site = 'Grass'
%site = 'JamesRes'
%site = 'LowDesert'
site = 'P301'
%site = 'Pinyon'
%site = 'PinyonBurn'
%site = 'Shorthair'
%site = 'SJER'
%site = 'Soaproot'
%site = 'Sage'

%from
yr1 = 2006;
%to
yr2 = 2012;
%------------------------------------------------------------------------%
%------------------------------------------------------------------------%
%ECsite
%------------------------------------------------------------------------%
ECroot = 'C:\towerData\WebsiteData\';
ECversion = '_v3_2.mat';

getthisECdata = [ECroot site ECversion];
load(getthisECdata);
%------------------------------------------------------------------------%
%set-up
%Cleaned-up best data - Use a handle
Dbest = DATA;

calm = Dbest(:,2) < 0.2;

%strike out the calm periods
Dbest(calm,20) = NaN;%Et
Dbest(calm,19) = NaN;%Fco2
Dbest(calm,18) = NaN;%LE
Dbest(calm,17) = NaN;%H

missingUst = isnan(Dbest(:,2)) == 1;

%strike out the missingUst periods
Dbest(missingUst,20) = NaN;%Et
Dbest(missingUst,19) = NaN;%Fco2
Dbest(missingUst,18) = NaN;%LE
Dbest(missingUst,17) = NaN;%H

%Fill the ET
%How many days to fill over?
if strcmp('JamesRes', site) == 1
    n_days_fill = 45;
else
   n_days_fill = 30; 
end

[Dbest, FilledET_regress, FilledET_robust, Kfill] = FillFh2o(Dbest, n_days_fill);


figure(1)
plot(Dbest(:,1), Kfill(:,1), 'k-')
ylabel('Kfill')

figure(2)
plot(Dbest(:,1), FilledET_regress(:,1), 'k-')
hold on
plot(Dbest(:,1), FilledET_robust(:,1), 'r-')
hold on
plot([2557 2557], [0 30], 'g-')
ylabel('ETfill')
legend('regress', 'robust');

%Fill the precip - assume no rain if data is missing
missing = isnan(DATA(:,8)) == 1;
DATA(missing,8) = 0;

%------------------------------------------------------------------------%
%find the year of the data
yr=365.2500;

y2006 = DATA(:,1) > 1         & DATA(:,1) <= yr +1;
y2007 = DATA(:,1) > yr +1     & DATA(:,1) <= (2*yr) +1;
y2008 = DATA(:,1) > (2*yr) +1 & DATA(:,1) <= (3*yr) +1;
y2009 = DATA(:,1) > (3*yr) +1 & DATA(:,1) <= (4*yr) +1;
y2010 = DATA(:,1) > (4*yr) +1 & DATA(:,1) <= (5*yr) +1;
y2011 = DATA(:,1) > (5*yr) +1 & DATA(:,1) <= (6*yr) +1;
y2012 = DATA(:,1) > (6*yr) +1 & DATA(:,1) <= (7*yr) +1;

%make a year colume
YR = DATA(:,1) * NaN;
YR(y2006,1) = 2006;
YR(y2007,1) = 2007;
YR(y2008,1) = 2008;
YR(y2009,1) = 2009;
YR(y2010,1) = 2010;
YR(y2011,1) = 2011;
YR(y2012,1) = 2012;
%------------------------------------------------------------------------%
%strike out 2010 at James because of bad sonic
if strcmp('JamesRes', site) == 1
    Dbest(y2010,20) = NaN;%Et
    Dbest(y2010,19) = NaN;%Fco2
    Dbest(y2010,18) = NaN;%LE
    Dbest(y2010,17) = NaN;%H
end

%------------------------------------------------------------------------%
%------------------------------------------------------------------------
%energy budget closure

%Is there soil_data
missing = sum(isnan(Dbest(:,16)));
lengthDbest16 = length(Dbest(:,16));

HMMMM_soildata = missing/lengthDbest16;


%availible energy
if HMMMM_soildata == 1
    AE= Dbest(:,15); %Rn - G
else
    AE= (Dbest(:,15) - Dbest(:,16)); %Rn - G
end

HLE = (Dbest(:,17) + Dbest(:,18)); %Sensible + Latent heat fluxes

isn = isnan(AE) | isnan(HLE);

ae=AE(~isn);
hle=HLE(~isn);
%------------------------------------------------------------------------%
%Here is Your energy budget closure forced through zero
energy_closure_zero = ae \ hle %Matlab finds the least square solution 
EB = (1/energy_closure_zero)
%------------------------------------------------------------------------%
%Here is the filled evaporation with correct units
Evap_regress = (FilledET_regress(:,1)*(1/1000)*18.01*(1/1000)*(1/1000)*1000*60*30) * EB; % mm/30min
Evap_robust = (FilledET_robust(:,1)*(1/1000)*18.01*(1/1000)*(1/1000)*1000*60*30) * EB; % mm/30min
%------------------------------------------------------------------------%
%sum up the annual ET
ET2006_regress = sum(Evap_regress(y2006,1));
ET2007_regress = sum(Evap_regress(y2007,1));
ET2008_regress = sum(Evap_regress(y2008,1));
ET2009_regress = sum(Evap_regress(y2009,1));
ET2010_regress = sum(Evap_regress(y2010,1));
ET2011_regress = sum(Evap_regress(y2011,1));
ET2012_regress = sum(Evap_regress(y2012,1));

%sum up the annual ET
ET2006_robust = sum(Evap_robust(y2006,1));
ET2007_robust = sum(Evap_robust(y2007,1));
ET2008_robust = sum(Evap_robust(y2008,1));
ET2009_robust = sum(Evap_robust(y2009,1));
ET2010_robust = sum(Evap_robust(y2010,1));
ET2011_robust = sum(Evap_robust(y2011,1));
ET2012_robust = sum(Evap_robust(y2012,1));

t2006 = mean(YR(y2006,1));
t2007 = mean(YR(y2007,1));
t2008 = mean(YR(y2008,1));
t2009 = mean(YR(y2009,1));
t2010 = mean(YR(y2010,1));
t2011 = mean(YR(y2011,1));
t2012 = mean(YR(y2012,1));

Time = [t2006; t2007; t2008; t2009; t2010; t2011; t2012];
ETs_regress=[ET2006_regress;ET2007_regress;ET2008_regress;ET2009_regress;ET2010_regress;ET2011_regress;ET2012_regress];
ETs_robust=[ET2006_robust;ET2007_robust;ET2008_robust;ET2009_robust;ET2010_robust;ET2011_robust;ET2012_robust];

ETsums = [Time, ETs_regress, ETs_robust];
   


