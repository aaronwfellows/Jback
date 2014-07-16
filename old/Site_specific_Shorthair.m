function [D, TEMP, H2O, rho_fix, v_wind] = Site_specific_Shorthair(HEADER, D, Press, R_mol, Tc, a)

%This file is a site specific clean-up for DMERGE data at the SJER in the Sierra Nevada Mountains
%awf 10/23/2012

%OUTPUTS:
%D is the filtered DMERGE data file
%TEMP is a cleaned and filled temperature in deg C 
%H2O is the water mixing ratio [mmol h2o/ mol of dry air]
%rho_fix is bad temp and water vapor measurements
%----> rho_fix is a handle that points to fast calculated
%----> fluxes that used bad dry density estimates 
%v_wind is cleaned up bad wind components

%IRGA calibrations
%The following is the data from the SiteDetails xls spreadsheet.  
%added the 30-Sep-2010 service date because there is a sharp change in the
%irga calibration - This will improve the h2o calibration
%IRGAchanged = [10,23,2009;9, 30, 2010; 9,29,2012]; 
IRGAchanged = [10,23,2009;9, 30, 2010];


if a==1

%most of the script is used to filter out bad data points
%generate a consensus temp, h2o, and wind vectors

%----------------------------------------------------------------------

%use ones through the program
len=length(D(1,:));
o=ones(1,len);

%use time to remove bad data
time = D(1,:);

%adjust the u wind dir
focus = time < 638;
D(207, focus) = D(207, focus)*100;

%additional pressure filtering
%focus = time < 605;
%bad = focus & D(110,:) < 55;
%D(110,bad) = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%By measurement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------
%HMP_Temp 
%----------------------------------------------------------------------
%data logger
badHMP_dl = isnan(D(97,:))==1 | isempty(D(97,:)) ==1 | isfinite(D(97,:)) == 0;
badHMP_dl = badHMP_dl | D(97,:)<-20;
badHMP_dl = badHMP_dl | time < 304.1;
D(97,badHMP_dl) = NaN;

%goes
badHMP_goes = isnan(D(223,:))==1 | isempty(D(223,:)) ==1 | isfinite(D(223,:)) == 0;
badHMP_goes = badHMP_goes | D(223,:)<-20;
badHMP_goes = badHMP_goes | time < 304.1;
D(223,badHMP_goes) = NaN;

%build a consensus HMP_Temp [C]
thmp=D(97,:);%HMP Temp
bad=isnan(thmp(1,:));
thmp(1,bad)=D(223,bad); %GOES TEMP

%----------------------------------------------------------------------
%Tsonic
%----------------------------------------------------------------------
%clean mat
Ts_mat_isnan=isnan(D(8,:))==1 | isempty(D(8,:)) ==1 | isfinite(D(8,:)) == 0;
bad_Ts_mat= Ts_mat_isnan ==1;
bad_Ts_mat= bad_Ts_mat | time < 304.13;
D(8,bad_Ts_mat)=NaN;

%clean dl
Ts_dl_isnan=isnan(D(84,:))==1 | isempty(D(84,:)) ==1 | isfinite(D(84,:)) == 0;
bad_Ts_dl = Ts_dl_isnan ==1;
bad_Ts_dl = bad_Ts_dl | time < 304.1;
D(84,bad_Ts_dl)=NaN;

%clean GOES
Ts_goes_isnan=isnan(D(210,:))==1 | isempty(D(210,:)) ==1 | isfinite(D(210,:)) == 0;
bad_Ts_goes = Ts_goes_isnan ==1;
focus = time > 342 & time < 343;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 402 & time < 403;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 475.5 & time < 477;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 483.15 & time < 484.8;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 715 & time < 722;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 777.5 & time < 778;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
bad_Ts_goes = bad_Ts_goes | time > 791.67 & time < 793;
focus = time > 795.87 & time < 797;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
bad_Ts_goes = bad_Ts_goes | time > 1038.2 & time < 1038.6;
D(210,bad_Ts_goes)=NaN;

%----------------------------------------------------------------------
%%   Water Vapor Pressure (HMP Probe) - convert RH to mixing ratio
%----------------------------------------------------------------------
%we cleaned up the hmp temp 
%==> fix up h2o_hmp calculation that were done in the data logger

%Lowe ~1977
A_0=6.107799961;
A_1=0.4436518521;
A_2=0.01428945805;
A_3=0.0002650648471;
A_4=0.000003031240396;
A_5=0.00000002034080948;
A_6=0.00000000006136820929;

%need to use the hmp temp because this is where the humidity is being
%measured
esat=0.1*(A_0+thmp.*(A_1+thmp.*(A_2+thmp.*(A_3+thmp.*(A_4+thmp.*(A_5+thmp.*A_6)))))); %[kPa]
 
%RH - clean-up relative humidity
%----------------------------------------------------------------------
%HMP_RH
%----------------------------------------------------------------------
%Work this out - awf 10/23/2012
%data logger
bad= D(98,:) < 0.01;
D(98,bad) = NaN;

%goes
bad= D(224,:) < 0.01;
bad= bad | D(224,:) > 20;
D(224,bad) = NaN;

e=D(98,:).*esat;%[kPa]
mol_frac_dl = e/Press; %mol h2o/mol wet air
mixing_dl=mol_frac_dl./(1-mol_frac_dl); %mol h2o/mol dry air
D(100,:)= mixing_dl*1000; %mmol h2o/mol dry air

e=D(224,:).*esat;%[kPa]
mol_frac_goes = e/Press;  
mixing_goes=mol_frac_goes./(1-mol_frac_goes); %mol h2o/mol dry air
D(226,:)=mixing_goes*1000; %mmol h2o/mol dry air

%----------------------------------------------------------------------
%IRGA H2O concentrations - mmol h20/mol dry air
%carry bad IRGA H2O concentrations  through the water fluxes
%----------------------------------------------------------------------
%clean fast
H2Oc_mat_isnan = isnan(D(42,:))==1 | isempty(D(42,:)) ==1 | isfinite(D(42,:)) == 0;
bad_H2Oc_mat = H2Oc_mat_isnan ==1;
bad_H2Oc_mat = bad_H2Oc_mat | D(42, :) > 29.8; %h2o conc are saturated at 30 mmol/mol 
D(42, bad_H2Oc_mat) = NaN;

%use this to remove the flux calculations
h2oflux_axe = bad_H2Oc_mat;

%there seems to be an offset problem not associated with calibrations
focus = time > 339 & time < 354;
focus = focus | time > 372.45 & time < 374.42;
focus = focus | time > 379 & time < 391.42;
focus = focus | time > 399.8 & time < 407;
focus = focus | time > 435 & time < 461;

focus_1=focus;
focus = focus & isfinite(D(42,:))==1;

p = polyfit(D(100,focus), D(42,focus),1); %p(1) is slope; p(2) is intercept

D(42, focus) = D(42, focus)/p(1);
%the Fh2o is also incorrect
D(47, focus) = D(47, focus)/p(1);
D(30, focus) = D(30, focus)/p(1);

%need to correct for bad water vapor concentration in the density
%calculation later
bad_H2Oc_mat= bad_H2Oc_mat | focus ==1;

%second half of the time series
focus = time > 679.7 & time < 685.82;
focus = focus | time > 695 & time < 696.38;
focus = focus | time > 700 & time < 702.38;
focus = focus | time > 710.4 & time < 711.4;
focus = focus | time > 716.5 & time < 717.4;
focus = focus | time > 763.4 & time < 764.38;
focus = focus | time > 783.4 & time < 784.39;
focus = focus | time > 790.5 & time < 791.5;
focus = focus | time > 793.6 & time < 794.38;
focus = focus | time > 797.6 & time < 798.38;
focus = focus | time > 881.4 & time < 882.4;
focus = focus | time > 971.36 & time < 973.4;
focus = focus | time > 1018.8 & time < 1019.38;
focus = focus | time > 1044.5 & time < 1045.38;
focus = focus | time > 1048.4& time < 1049.4;
focus = focus | time > 1099.4 & time < 1100.4;
focus = focus | time > 1105.86 & time < 1106.4;
focus = focus | time > 1120.45 & time < 1121.4;

focus_2=focus;
focus = focus & isfinite(D(42,:))==1;

p = polyfit(D(100,focus), D(42,focus),1); %p(1) is slope; p(2) is intercept

%need to correct for bad water vapor concentration in the density calculation later
bad_H2Oc_mat= bad_H2Oc_mat | focus ==1;

D(42, focus) = D(42, focus)/p(1);
%the Fh2o is also incorrect
D(47, focus) = D(47, focus)/p(1);
D(30, focus) = D(30, focus)/p(1);

%clean data logger
H2Oc_dl_isnan = isnan(D(86,:))==1 | isempty(D(86,:)) ==1 | isfinite(D(86,:)) == 0;
bad_H2Oc_dl = H2Oc_dl_isnan ==1;
bad_H2Oc_dl = bad_H2Oc_dl | D(86, :) > 29.8; %h2o conc are saturated at 30 mmol/mol 
bad_H2Oc_dl = bad_H2Oc_dl | D(78, :) < 0; %no variabiltiy 

D(86, bad_H2Oc_dl) = NaN;

focus = focus_1 & isfinite(D(100,:)) & isfinite(D(86,:));
p = polyfit(D(100,focus), D(86,focus),1); %p(1) is slope; p(2) is intercept
D(86, focus) = D(86, focus)/p(1);
%may have to convert uz_h2o for Fh2o calculation
D(64, focus) = D(64, focus)/p(1);


focus = focus_2 & isfinite(D(100,:)) & isfinite(D(86,:));
p = polyfit(D(100,focus), D(86,focus),1); %p(1) is slope; p(2) is intercept
D(86, focus) = D(86, focus)/p(1);
%may have to convert uz_h2o for Fh2o calculation
D(64, focus) = D(64, focus)/p(1);

%convert to mixing ratio
mixing_dl=D(86, :)./(1-(D(86, :)/1000)); %mol h2o/mol dry air
D(86, :)= mixing_dl; %mmol h2o/mol dry air

%clean GOES
H2Oc_goes_isnan = isnan(D(212,:))==1 | isempty(D(212,:)) ==1 | isfinite(D(212,:)) == 0;
bad_H2Oc_goes = H2Oc_goes_isnan ==1;
bad_H2Oc_goes = bad_H2Oc_goes | D(212, :) > 29.8; %h2o conc are saturated at 30 mmol/mol 
bad_H2Oc_goes = bad_H2Oc_goes | D(204, :) == 0; %no variabilty 
bad_H2Oc_goes = bad_H2Oc_goes | D(212, :) < -200; %no variabilty 
D(212, bad_H2Oc_goes) = NaN;

focus = focus_1 & isfinite(D(100,:)) & isfinite(D(212,:));
p = polyfit(D(100,focus), D(212,focus),1); %p(1) is slope; p(2) is intercept
D(212, focus) = D(212, focus)/p(1);
D(190, focus) = D(190, focus)/p(1);

focus = focus_2 & isfinite(D(100,:)) & isfinite(D(212,:));
p = polyfit(D(100,focus), D(212,focus),1); %p(1) is slope; p(2) is intercept
D(212, focus) = D(212, focus)/p(1);
%may have to convert uz_h2o for Fh2o calculation
D(190, focus) = D(190, focus)/p(1);

%convert to mixing ratio
mixing_goes=D(212, :)./(1-(D(212, :)/1000)); %mol h2o/mol dry air
D(212, :)= mixing_goes; %mmol h2o/mol dry air

%----------------------------------------------------------------------
%build a consensus H2O mixing ratio [mmol h20/mol of dry air]
H2O=D(100,:);%HMP 
bad=isnan(H2O(1,:));
H2O(1,bad)=D(226,bad); %GOES HMP
bad=isnan(H2O(1,:));
H2O(1,bad)=D(42,bad);%fast
bad=isnan(H2O(1,:));
H2O(1,bad)=D(86,bad);%dl
bad=isnan(H2O(1,:));
H2O(1,bad)=D(212,bad);%goes

%a fill H2O - h2o has a small effect on density ~ few percent - so we
%should fill in h2o for density corrections 
medianh2o=nanmedian(H2O(1,:));
h2ofill=medianh2o.*o;

bad=isnan(H2O(1,:));
H2O(1,bad)=h2ofill(1,bad);

%----------------------------------------------------------------------
%convert sonic temp over to actual temperature and fill
%----------------------------------------------------------------------
%water mixing ratio to mol fraction [mmol h20/mol of wet air]
h2o_mol_frac = H2O./(1+(H2O/1000));
%----------------------------------------------------------------------
%Kaimal and Gaynor - 1991
%Tkact_mat   = (D(8,:))  ./ (1 + 0.00032 * h2o_mol_frac); %[C]
%Tkact_dl   = (D(84,:))  ./ (1 + 0.00032 * h2o_mol_frac); %[C]
%Tkact_goes   = (D(210,:))  ./ (1 + 0.00032 * h2o_mol_frac); %[C]

Tkact_mat   = (D(8,:)+273.15)  ./ (1 + 0.00032 * h2o_mol_frac); %[K]
Tkact_dl   = (D(84,:)+273.15)  ./ (1 + 0.00032 * h2o_mol_frac); %[K]
Tkact_goes   = (D(210,:)+273.15)  ./ (1 + 0.00032 * h2o_mol_frac); %[K]

Tkact_mat   = Tkact_mat-273.15; %[C]
Tkact_dl   = Tkact_dl-273.15; %[C]
Tkact_goes   = Tkact_goes-273.15; %[C]
%----------------------------------------------------------------------
%build a consensus Temp
TEMP=D(97,:);%HMP Temp
bad=isnan(TEMP(1,:));
TEMP(1,bad)=D(223,bad); %GOES TEMP
bad=isnan(TEMP(1,:));
TEMP(1,bad)=Tkact_mat(1,bad);%fast
bad=isnan(TEMP(1,:));
TEMP(1,bad)=Tkact_dl(1,bad);%dl
bad=isnan(TEMP(1,:));
TEMP(1,bad)=Tkact_goes(1,bad);%goes
%----------------------------------------------------------------------

%fast calculated a bad dry density if the water vapor was wrong or temp was
%wrong - The error was probably small but we can adjust this with the ratio
%of best_dry_density/dry_density_that_was_used ==> Rho Calc in combined ./ D(48,:)

%When rho_fix ==1, we should apply the density correction
rho_fix = bad_Ts_mat | bad_H2Oc_mat;
%~49% of the observations will be adjusted
%----------------------------------------------------------------------
%winds
%----------------------------------------------------------------------
%verticle windspeed (w) is incorrect for all time < 740 for all data sources

%component wind problems - variance is 0 is impossible
%u_wind
u_mat_var=D(9,:);
Umat_isnan=isnan(D(5,:))==1 | isempty(D(5,:)) ==1 | isfinite(D(5,:)) == 0;
ubad_wind_mat=u_mat_var(1,:)==0;
ubad_wind_mat = ubad_wind_mat | Umat_isnan(1,:)==1;
ubad_wind_mat = ubad_wind_mat | time > 1177.59 & time < 1233; % tower tilted
D(5,ubad_wind_mat)=NaN;

%v_wind
v_mat_var=D(10,:);
Vmat_isnan=isnan(D(6,:))==1 | isempty(D(6,:)) ==1 | isfinite(D(6,:)) == 0;
vbad_wind_mat = v_mat_var(1,:)==0;
vbad_wind_mat = vbad_wind_mat | Vmat_isnan(1,:)==1;
vbad_wind_mat = vbad_wind_mat | time > 1177.59 & time < 1233;
D(6,vbad_wind_mat)=NaN;

w_mat_var=D(11,:);
Wmat_isnan=isnan(D(7,:))==1 | isempty(D(7,:)) ==1 | isfinite(D(7,:)) == 0;
wbad_wind_mat = w_mat_var(1,:)==0;
wbad_wind_mat = wbad_wind_mat | Wmat_isnan(1,:)==1;
wbad_wind_mat = wbad_wind_mat | time > 1177.59 & time < 1233;
D(7,wbad_wind_mat)=NaN;


%clean dl
%DATA LOGGER
u_dl_var=D(66,:);%'Ux variance
Udl_isnan=isnan(D(81,:))==1 | isempty(D(81,:)) ==1 | isfinite(D(81,:)) == 0;
ubad_wind_dl=u_dl_var(1,:)==0;
ubad_wind_dl = ubad_wind_dl | Udl_isnan(1,:)==1;
ubad_wind_dl = ubad_wind_dl | time > 1177.59 & time < 1233;
D(81,ubad_wind_dl)=NaN;

v_dl_var=D(71,:);%'Uy variance
Vdl_isnan=isnan(D(82,:))==1 | isempty(D(82,:)) ==1 | isfinite(D(82,:)) == 0;
vbad_wind_dl = v_dl_var(1,:)==0;
vbad_wind_dl = vbad_wind_dl | Vdl_isnan(1,:)==1;
vbad_wind_dl = vbad_wind_dl | time > 1177.59 & time < 1233;
D(82,vbad_wind_dl)=NaN;

w_dl_var=D(60,:);% Uz variance
Wdl_isnan=isnan(D(83,:))==1 | isempty(D(83,:)) ==1 | isfinite(D(83,:)) == 0;
wbad_wind_dl = w_dl_var(1,:)==0;
wbad_wind_dl = wbad_wind_dl | Wdl_isnan(1,:)==1;
wbad_wind_dl = wbad_wind_dl | time > 1177.59 & time < 1233;
D(83,wbad_wind_dl)=NaN;


%clean goes
u_goes_var = D(192,:);% Ux variance
Ugoes_isnan=isnan(D(207,:))==1 | isempty(D(207,:)) ==1 | isfinite(D(207,:)) == 0;
ubad_wind_goes = u_goes_var(1,:) == 0;
ubad_wind_goes = ubad_wind_goes | Ugoes_isnan ==1;
ubad_wind_goes = ubad_wind_goes | time > 1177.59 & time < 1233;
D(207,ubad_wind_goes) = NaN;% Ux

v_goes_var = D(197,:);% Uy variance 
Vgoes_isnan=isnan(D(208,:))==1 | isempty(D(208,:)) ==1 | isfinite(D(208,:)) == 0;
vbad_wind_goes = v_goes_var==0;
vbad_wind_goes = vbad_wind_goes | Vgoes_isnan ==1;
vbad_wind_goes = vbad_wind_goes | time > 1177.59 & time < 1233;
D(208,vbad_wind_goes) = NaN;% Uy

w_goes_var = D(186,:);% Uz variance
Wgoes_isnan=isnan(D(209,:))==1 | isempty(D(209,:)) ==1 | isfinite(D(209,:)) == 0;
wbad_wind_goes = w_goes_var==0;
wbad_wind_goes = wbad_wind_goes | Wgoes_isnan ==1;
wbad_wind_goes = wbad_wind_goes | time > 1177.59 & time < 1233;
D(209,wbad_wind_goes) = NaN;% Uy

%fill winds
%u wind
u=D(5,:);%u wind
bad=isnan(u(1,:));
u(1,bad)=D(81,bad); %u wind
bad=isnan(u(1,:));
u(1,bad)=D(207,bad);%u wind

%v wind
v=D(6,:);%v wind
bad=isnan(v(1,:));
v(1,bad)=D(82,bad); %v wind
bad=isnan(v(1,:));
v(1,bad)=D(208,bad);%v wind

%w wind
w=D(7,:);%w wind
bad=isnan(w(1,:));
w(1,bad)=D(83,bad); %w wind
bad=isnan(w(1,:));
w(1,bad)=D(209,bad);%w wind

v_wind=[u;v;w];

%----------------------------------------------------------------------
%Wind Direction
%----------------------------------------------------------------------
%mat wind dir is wrong when component winds are wrong
bad_wind_mat = ubad_wind_mat | vbad_wind_mat | wbad_wind_mat;
D(22,bad_wind_mat)=NaN;

%----------------------------------------------------------------------
%Ustar
%----------------------------------------------------------------------
%clean FAST
bad = D(9,:) == 0; % u variance == 0 
bad = bad | isnan(D(9,:)) == 1; %missing u
bad = bad | D(10,:) == 0; % u variance == 0 
bad = bad | isnan(D(10,:)) == 1; %missing u
bad = bad | D(11,:) == 0; %uz variability = 0
bad = bad | isnan(D(11,:)) == 1; %missing uz
bad = bad | time > 1177.59 & time < 1233;
D(21,bad)=NaN;

%clean data logger
%bad = bad_wind_dl;
bad = D(60,:) == 0;
bad = bad | isnan(D(60,:)) == 1; % is there uz uz data
bad = bad | D(66,:) == 0; % sensor not working
bad = bad | isnan(D(66,:)) == 1; % is there ux ux data
bad = bad | time > 1177.59 & time < 1233;
D(61,bad)=NaN;

bad = D(60,:) == 0;
bad = bad | isnan(D(60,:)) == 1; % is there uz uz data
bad = bad | D(71,:) == 0; % sensor not working
bad = bad | isnan(D(71,:)) == 1; % is there uy uy data
bad = bad | time > 1177.59 & time < 1233;
D(62,bad)=NaN;

%clean GOES
bad = D(186,:) == 0; % uz uz
bad = bad | isnan(D(186,:)) == 1; % is there uz uz data
bad = bad | D(192,:) == 0; % sensor not working
bad = bad | isnan(D(192,:)) == 1; % is there ux ux data
bad = bad | time > 1177.59 & time < 1233;
D(187,bad)=NaN;

bad = D(186,:) == 0; % uz uz
bad = bad | isnan(D(186,:)) == 1; % is there uz uz data
bad = bad | D(197,:) == 0; % sensor not working
bad = bad | isnan(D(197,:)) == 1; % is there uy uy data
bad = bad | time > 1177.59 & time < 1233;
D(188,bad)=NaN;

%----------------------------------------------------------------------
%Sensible heat
%----------------------------------------------------------------------
%clean fast - this is despiked in fastflux
bad = D(12,:) == 0; % Ts var == 0 
bad = bad | isnan(D(12,:)) == 1; %missing Ts
bad = bad | D(11,:) == 0; %w variability = 0
bad = bad | isnan(D(11,:)) == 1; %missing w
bad = bad | time > 1177.59 & time < 1233;
D(29,bad)=NaN;

bad = D(29,:) < -300;
bad = bad | D(29,:) > 800;
D(29,bad)=NaN;

%clean data logger
%need a bad covariance term
bad = D(60,:) == 0;
bad = bad | isnan(D(60,:)) == 1; % is there uz uz data
bad = bad | isnan(D(80,:)) == 1; % is there uy uy data
bad = bad | D(80,:) == 0; % sensor not working
bad = bad | time > 1177.59 & time < 1233;
D(65,bad)=NaN;


%clean GOES 
bad = D(186,:) == 0; % uz uz
bad = bad | isnan(D(186,:)) == 1; % is there uz uz data
bad = bad | isnan(D(206,:)) == 1; % is there ux ux data
bad = bad | D(206,:) == 0; % sensor not working
bad = bad | time > 1177.59 & time < 1233;
D(191, bad) = NaN;

%----------------------------------------------------------------------
%CO2 concentration
%----------------------------------------------------------------------
%clean fast
bad_fastCO2 = time > 1177.59 & time < 1233; 
D(35,bad_fastCO2)=NaN;

%clean data logger
bad_dlCO2 = time > 1177.59 & time < 1233;
D(85, bad_dlCO2) = NaN;

D(85,:)=D(85,:)./(1-(h2o_mol_frac(1,:)/1000)); %co2/mol dry air

%an additional filter on co2 conc - does not remove the Fco2
bad = D(135,:) > 0.01;
D(85, bad) = NaN;

%clean GOES
bad_goesCO2 = time > 1177.59 & time < 1233;
D(211, bad_goesCO2) = NaN;

D(211,:)=D(211,:)./(1-(h2o_mol_frac(1,:)/1000)); %co2/mol dry air

%an additional filter on co2 conc - does not remove the Fco2
bad = D(261,:) > 0.01;
D(211, bad) = NaN;
%----------------------------------------------------------------------
%FCO2
%----------------------------------------------------------------------
%clean fast
bad = bad_fastCO2;
bad = bad | D(11,:) == 0; %w variability = 0
bad = bad | isnan(D(11,:)) == 1; %missing w
bad = bad | D(36,:) == 0; %co2 variability = 0
bad = bad | isnan(D(36,:)) == 1; %missing co2
D(46, bad) = NaN;

%out of reasonable range range
bad = D(46,:) < -20;
bad = bad | D(46,:) > 10;
D(46, bad) = NaN;

%clean data logger FCO2 by taking out uz_co2 covariance
bad_fco2_dl = bad_dlCO2;
bad_fco2_dl = bad_fco2_dl | isnan(D(60,:))==1; % no uz uz
bad_fco2_dl = bad_fco2_dl | D(60,:)==0;%no wind variability
bad_fco2_dl = bad_fco2_dl | D(75,:)==0; % co2 variability
bad_fco2_dl = bad_fco2_dl | isnan(D(75,:)) == 1; % no data
bad_fco2_dl = bad_fco2_dl | isnan(D(63,:)) == 1; % no data
D(63, bad_fco2_dl) = NaN;


%clean GOES
%bad_fco2_goes = bad_goesCO2;
bad_fco2_goes = isnan(D(186,:))==1; % sonic is not working no uz uz
bad_fco2_goes = bad_fco2_goes | D(186,:)==0;%no wind variability
bad_fco2_goes = bad_fco2_goes | D(201,:)==0; % co2 variability
bad_fco2_goes = bad_fco2_goes | isnan(D(201,:))==1; % sonic is not working no uz uz
bad_fco2_goes = bad_fco2_goes | isnan(D(189,:)) == 1; % no data
D(189, bad_fco2_goes) = NaN;

%----------------------------------------------------------------------
%Adjust water fluxes using IRGA H2O variabilty compared with HMP variability
%NOTE: IRGA H2O variabilty seems high compared with HMP variability.
%Therefore, we over estimate FH2O.  This will reduce the water flux.

%Make sure the water vapor concentrations are cleaned first.
%----------------------------------------------------------------------
Cnums=datenum(IRGAchanged(:,3), IRGAchanged(:,1), IRGAchanged(:,2));
Cnumsince2009=Cnums-datenum(2009,1,0); % this is for the Shorthair

firstday = min(time);
lastday = max(time);

length_Cnumsince2009=length(Cnumsince2009(:,1));
Cnumsince2009(1,1)=firstday-1;
Cnumsince2009(length_Cnumsince2009+1,1)=lastday;
length_Cnumsince2009 = length_Cnumsince2009+1;

% serviced figure
f1=(Cnumsince2009./Cnumsince2009)*40;
f2=(Cnumsince2009./Cnumsince2009)*-5;
f3=(Cnumsince2009./Cnumsince2009)*40;
F=[f1,f2,f3]';
[m,n]=size(F);
l=m*n;
F=reshape(F,l,1);
C09=[Cnumsince2009; Cnumsince2009; Cnumsince2009];
C09=sort(C09);

figure(2)
plot(time, D(100,:), 'k.');
hold on
plot(time, D(42,:), 'b.');
hold on
plot(time, D(126,:), 'g.');
hold on
plot(C09, F, 'r-');
legend('black=hmp', 'blue=irga', 'green=h2o zero (using span gas', 'irga calibration');

Pmat=[];
Pdl=[];
Pgoes=[];


for i=2:length_Cnumsince2009
    %proc data
    good_mat = time(1,:) < Cnumsince2009(i)-2 & time(1,:) > Cnumsince2009(i-1)+2 & isfinite(D(42,:))==1 & isfinite(D(100,:))==1;
    irga=D(42,good_mat);
    hmp=D(100,good_mat);
    p = polyfit(hmp,irga,1) %p(1) is slope; p(2) is intercept
    c_tm = time(1,:) <= Cnumsince2009(i) & time(1,:) > Cnumsince2009(i-1); %adjustment period
    
    %record so you can look later
    tt= median(time(1,good_mat)) 
    Pmat=[Pmat; tt, p];
    
    %The offset does not make sense here. Adjust water fluxes by dividing by the slope.
    outrageous_p = isnan(p(1,1))==1 | isfinite(p(1,1)) == 0 | p(1,1) > 1.2 | p(1,1) < 0.8 | length(irga) < 30
    
    figure(1)
    plot(hmp,irga, 'k.')
    hold on
    plot(hmp, p(1,1)*hmp+p(1,2), 'r-')
    xlabel('hmp h2o')
    ylabel('mat irga h2o')
    pause 
    close(figure(1))
    
    if outrageous_p == 1
        disp('can not adjust h2o fluxes')
        tt
    else
        D(47,c_tm)=D(47,c_tm)./p(1,1); %Fh20
        D(30,c_tm)=D(30,c_tm)./p(1,1); %LH
    end
      
    
    %data logger
    good_dl = time(1,:) < Cnumsince2009(i)-2 & time(1,:) > Cnumsince2009(i-1)+2 & isfinite(D(86,:))==1 & isfinite(D(100,:))==1;
    irga=D(86,good_dl);
    hmp=D(100,good_dl);
    p = polyfit(hmp,irga,1) %p(1) is slope; p(2) is intercept
    c_tm = time(1,:) <= Cnumsince2009(i) & time(1,:) > Cnumsince2009(i-1); %adjustment period
    
    %record so you can look later
    tt= median(time(1,good_dl)) 
    Pdl=[Pdl; tt, p];
    
    outrageous_p = isnan(p(1,1))==1 | isfinite(p(1,1)) == 0 | p(1,1) > 1.2 | p(1,1) < 0.8 | length(irga) < 30
    
    figure(1)
    plot(hmp,irga, 'k.')
    hold on
    plot(hmp, p(1,1)*hmp+p(1,2), 'r-')
    xlabel('hmp h2o')
    ylabel('dl irga h2o')
    pause 
    close(figure(1))
    
    if outrageous_p == 1
        disp('can not adjust h2o fluxes')
        tt
    else
    %The offset does not make sense here. Adjust water fluxes by dividing by the slope.
        D(64,c_tm)=D(64,c_tm)./p(1,1); %cov h2o dl
    end
    
    
    %goes
    good_goes = time(1,:) < Cnumsince2009(i)-2 & time(1,:) > Cnumsince2009(i-1)+2 & isfinite(D(212,:))==1 & isfinite(D(226,:))==1;
    irga=D(212,good_goes);
    hmp=D(226,good_goes);
    p = polyfit(hmp,irga,1) %p(1) is slope; p(2) is intercept
    c_tm = time(1,:) <= Cnumsince2009(i) & time(1,:) > Cnumsince2009(i-1); %adjustment period
    
    %record so you can look later
    tt= median(time(1,good_goes)) 
    Pgoes=[Pgoes; tt, p];
    
    %The offset does not make sense here. Adjust water fluxes by dividing by the slope.
    outrageous_p = isnan(p(1,1))==1 | isfinite(p(1,1)) == 0 | p(1,1) > 1.2 | p(1,1) < 0.8 | length(irga) < 30
    
    figure(1)
    plot(hmp,irga, 'k.')
    hold on
    plot(hmp, p(1,1)*hmp+p(1,2), 'r-')
    xlabel('hmp h2o')
    ylabel('goes irga h2o')
    pause 
    close(figure(1))
    
    if outrageous_p == 1
        disp('can not adjust h2o fluxes')
        tt
    else
        D(190,c_tm)=D(190,c_tm)./p(1,1); %cov h2o goes
    
    end
    
end

close(figure(2))
save('C:\towerData\combined\CombinedInfo\Shorthair\HMP_h2o_calibration_Shorthair.mat', 'Pmat', 'Pdl', 'Pgoes');

%----------------------------------------------------------------------
%FH2O
%----------------------------------------------------------------------
%clean fast

bad_mat_et = h2oflux_axe; %
bad_mat_et = bad_mat_et | D(43,:) == 0; %
bad_mat_et = bad_mat_et | isnan(D(43,:)) == 1; %
bad_mat_et = bad_mat_et | D(11,:) == 0; %uz variability = 0
bad_mat_et = bad_mat_et | isnan(D(11,:)) == 1; %missing uz
bad_mat_et = bad_mat_et | D(47,:) < -10;
bad_mat_et = bad_mat_et | D(47,:) > 15;
bad_mat_et = bad_mat_et | time > 1177.59 & time < 1233;

D(47, bad_mat_et) = NaN;

%clean data logger
bad_dl_et = bad_H2Oc_dl;
bad_dl_et = bad_dl_et | isnan(D(60,:))==1; % sonic is not working no uz uz
bad_dl_et = bad_dl_et | D(60,:)==0;%no wind variability
bad_dl_et = bad_dl_et | D(78,:)==0; % h2o variability
bad_dl_et = bad_dl_et | isnan(D(78,:))==1; % sonic is not working no uz uz
bad_dl_et = bad_dl_et | isnan(D(64,:)) == 1; % no data
bad_dl_et = bad_dl_et | time > 1177.59 & time < 1233;

D(64, bad_dl_et) = NaN;

%clean GOES
bad_goes_et = bad_H2Oc_goes;
bad_goes_et = bad_goes_et | isnan(D(186,:))==1; % sonic is not working no uz uz
bad_goes_et = bad_goes_et | D(186,:)==0;%no wind variability
bad_goes_et = bad_goes_et | D(204,:)==0; % h2o variability
bad_goes_et = bad_goes_et | isnan(D(204,:))==1; % 
bad_goes_et = bad_goes_et | isnan(D(190,:)) == 1; % no data
bad_goes_et = bad_goes_et | time > 1177.59 & time < 1233;

D(190, bad_goes_et) = NaN;

%----------------------------------------------------------------------
%Latent heat
%----------------------------------------------------------------------
%fast
%bad_le_mat = bad_mat_et | D(30,:) == 0;
bad_le_mat = bad_mat_et;
D(30,bad_le_mat)=NaN;

%----------------------------------------------------------------------
%Rn
%----------------------------------------------------------------------
%clean data logger
bad = time > 1177 & time < 1235; 
D(87,bad) = NaN;

%clean goes
bad = time > 1177 & time < 1235; 
D(213,bad) = NaN;

%----------------------------------------------------------------------
%Incoming solar radiation - pyranometer = SOLAR_IN
%----------------------------------------------------------------------
%clean data logger
bad = time > 1177 & time < 1235; 
D(88,bad) = NaN;

%clean goes
bad = time > 1177 & time < 1235; 
D(214,bad) = NaN;

%----------------------------------------------------------------------
%Outgoing solar radiation - pyranometer = SOLAR_OUT
%----------------------------------------------------------------------
%clean data logger
bad = time > 1177 & time < 1235; 
D(89,bad) = NaN;

%clean goes
bad = time > 1177 & time < 1235; 
D(215,bad) = NaN;

%----------------------------------------------------------------------
%PAR_In
%----------------------------------------------------------------------
%clean data logger
bad = time > 1177 & time < 1235; 
D(90,bad) = NaN;

%clean goes
bad = time > 1177 & time < 1235; 
D(216,bad) = NaN;

%----------------------------------------------------------------------
%PAR_Out
%----------------------------------------------------------------------
%clean data logger
bad = time > 1177 & time < 1235; 
D(91,bad) = NaN;

%clean goes
bad = time > 1177 & time < 1235; 
D(217,bad) = NaN;

%----------------------------------------------------------------------
%T_107
%----------------------------------------------------------------------
%dl
bad = D(102,:) < -30;
D(102,bad) = NaN;

%goes
bad = D(228,:) < -30;
bad = bad | time < 638.82;
D(228,bad) = NaN;

%----------------------------------------------------------------------
%Rain
%----------------------------------------------------------------------
%data logger
bad=D(103,:)>5000;
D(103,bad) = NaN;

%goes
bad=D(229,:)>200;
D(229,bad) = NaN;
%----------------------------------------------------------------------
%NDVI
%----------------------------------------------------------------------
bad = D(96,:) > 1;
D(96,bad)=NaN;

bad = D(222,:) > 1;
D(222,bad)=NaN;

%----------------------------------------------------------------------
%Soil T
%----------------------------------------------------------------------

%----------------------------------------------------------------------
%Clean Soil Moisture
%----------------------------------------------------------------------

%----------------------------------------------------------------------
%Fuel_M
%----------------------------------------------------------------------


%----------------------------------------------------------------------
%LWS
%----------------------------------------------------------------------


%----------------------------------------------------------------------
%Clean Matric potential sensors 
%----------------------------------------------------------------------


elseif a==2
  
%Mop-up

TEMP =[];
H2O =[];
rho_fix =[];
v_wind = [];   
  
%Sensible Heat
%dl
bad = D(339,:)<-300;
bad = bad | D(339,:)> 800;
D(339,bad)=NaN;

%goes
bad = D(340,:)<-300;
bad = bad | D(340,:)> 800;
D(340, bad) = NaN;

%Fh2o
%dl
bad_dl_et = D(347,:) < -10;
bad_dl_et = bad_dl_et | D(347,:) > 15;
D(347, bad_dl_et) = NaN;

%goes
bad_goes_et = D(348,:) < -10;
bad_goes_et = bad_goes_et | D(348,:) > 15;
D(348, bad_goes_et) = NaN;

%latent heat
%clean data logger
%bad_le_dl = bad_dl_et | D(343,:) == 0;
bad_le_dl = bad_dl_et;
D(343,bad_le_dl)=NaN;

%clean GOES 
%bad_le_goes = bad_goes_et | D(344,:) == 0;
bad_le_goes = bad_goes_et;
D(344,bad_le_goes)=NaN;

%FCO2
%dl
bad = D(345,:) < -20;
bad = bad | D(345,:) > 10;
D(345, bad) = NaN;

%goes
bad = D(346,:) < -20;
bad = bad | D(346,:) > 10;
D(346, bad) = NaN;

%Ustar
%clean dl
bad = D(322,:) > 2.5;
D(322,bad)=NaN;

%clean GOES
bad = D(323,:) > 2.5;
D(323,bad)=NaN;

%There are a few outrageous windspeeds 
bad = D(328,:) > 12;
D(328,bad)=NaN;

bad = D(329,:) > 12;
D(329,bad)=NaN;
end