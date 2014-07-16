function [D, TEMP, H2O, rho_fix, v_wind] = Site_specific_DC_LowDesWorking(HEADER, D, Press, R_mol, Tc, a)

%When were the calibration changes for the irga?
IRGAchanged = [4,21,2006;1,17,2007;2,28,2008;3,30,2010;5,4,2010;8,25,2011;6,22,2012];

%This file is a site specific clean-up for DMERGE data at the Low Desert
%awf 10/19/2012

%OUTPUTS:
%D is the filtered DMERGE data file
%TEMP is a cleaned and filled temperature in deg C 
%H2O is the water mixing ratio [mmol h2o/ mol of dry air]
%rho_fix is bad temp and water vapor measurements
%----> rho_fix is a handle that points to fast calculated
%----> fluxes that used bad dry density estimates 
%v_wind is cleaned up bad wind components


if a==1
%----------------------------------------------------------------------
%fix up the goes wind problem in the u direction
focus= D(1,:) > 502.22 & D(1,:) < 516.4;
focus= focus | D(1,:) > 529.26 & D(1,:) < 538.2;
focus= focus | D(1,:) > 1288.13;
D(207,~focus)= D(207,~focus).*100; %this is divided 
%----------------------------------------------------------------------

%use ones through the program
len=length(D(1,:));
o=ones(1,len);

%use time to remove bad data
time = D(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time stamp problem -- kill the goes data until this is figured out
bad = time > 2593 & time < 2637;
D(181:275,bad)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%By measurement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------
%HMP_Temp 
%----------------------------------------------------------------------
%data logger
badHMP_dl = isnan(D(97,:))==1 | isempty(D(97,:)) ==1 | isfinite(D(97,:)) == 0;
badHMP_dl = badHMP_dl | D(97,:) < -15;
badHMP_dl = badHMP_dl | time > 287.43 & time < 300.975; 
badHMP_dl = badHMP_dl | time > 303.34 & time < 312.93; 
badHMP_dl = badHMP_dl | time > 678.43 & time < 789.79; 
badHMP_dl = badHMP_dl | time > 887 & time < 1075.025; 
badHMP_dl = badHMP_dl | time > 1088.7 & time < 1090.99; 
badHMP_dl = badHMP_dl | time > 1114 & time < 1159.87; 
badHMP_dl = badHMP_dl | time > 1171.78 & time < 1288; 
badHMP_dl = badHMP_dl | time > 1483.42 & time < 1484.436;
D(97,badHMP_dl) = NaN;

%goes
badHMP_goes = isnan(D(223,:))==1 | isempty(D(223,:)) ==1 | isfinite(D(223,:)) == 0;
badHMP_goes = badHMP_goes | D(223,:) < -15;
badHMP_goes = badHMP_goes | time > 118.81 & time < 118.82; 
badHMP_goes = badHMP_goes | time > 287.43 & time < 300.975; 
badHMP_goes = badHMP_goes | time > 303.34 & time < 312.93;
badHMP_goes = badHMP_goes | time > 678.43 & time < 789.79; 
badHMP_goes = badHMP_goes | time > 836.83 & time < 836.86; 
badHMP_goes = badHMP_goes | time > 887 & time < 1075.025; 
badHMP_goes = badHMP_goes | time > 1088.7 & time < 1090.99; 
badHMP_goes = badHMP_goes | time > 1114 & time < 1159.87; 
badHMP_goes = badHMP_goes | time > 1171.78 & time < 1288.1; 
badHMP_goes = badHMP_goes | time > 1483.42 & time < 1484.436;
%badHMP_goes = badHMP_goes | time > 2592 & time < 2638;
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
bad_Ts_mat = bad_Ts_mat | time > 332.74 & time < 332.76;
bad_Ts_mat = bad_Ts_mat | time > 382.97 & time < 382.99;
bad_Ts_mat = bad_Ts_mat | time > 1861.05 & time < 1861.07;
D(8,bad_Ts_mat)=NaN;

%clean dl
Ts_dl_isnan=isnan(D(84,:))==1 | isempty(D(84,:)) ==1 | isfinite(D(84,:)) == 0;
bad_Ts_dl= Ts_dl_isnan ==1;
D(84,bad_Ts_dl)=NaN;

%clean GOES
Ts_goes_isnan=isnan(D(210,:))==1 | isempty(D(210,:)) ==1 | isfinite(D(210,:)) == 0;
bad_Ts_goes = Ts_goes_isnan ==1 | D(210,:)==0; 
bad_Ts_goes = bad_Ts_goes | D(210,:)< -100; 
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
 
%RH 
bad = time > 861.96 & time < 871.41;
bad = bad | time > 1171.68 & time < 1228;
bad = bad | time > 1273.78 & time < 1275;
focus = time > 1483 & time < 1485;
bad = bad | focus & D(98,:) < 0.01;
D(98,bad)=NaN;

bad = time > 861.96 & time < 871.41;
bad = bad | D(224,:) < -50;
bad = bad | time > 207.854 & time < 207.8545;
bad = bad | time > 1171.68 & time < 1228;
bad = bad | time > 1273.78 & time < 1275;
focus = time > 1062 & time < 1064;
bad = bad | focus & D(224,:)==0;
focus = time > 700 & time < 701;
bad = bad | focus & D(224,:)==0;
focus = time > 1483 & time < 1485;
bad = bad | focus & D(224,:)==0;
D(224,bad)=NaN;

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
bad_H2Oc_mat = bad_H2Oc_mat | D(42, :) < 0.2; %h2o conc are saturated at 30 mmol/mol 
bad_H2Oc_mat = bad_H2Oc_mat | D(43, :) == 0; % no variabilty
D(42, bad_H2Oc_mat) = NaN;

%clean data logger
H2Oc_dl_isnan = isnan(D(86,:))==1 | isempty(D(86,:)) ==1 | isfinite(D(86,:)) == 0;
bad_H2Oc_dl = H2Oc_dl_isnan ==1;
bad_H2Oc_dl = bad_H2Oc_dl | D(86, :) > 29.8; %h2o conc are saturated at 30 mmol/mol 
bad_H2Oc_dl = bad_H2Oc_dl | D(86, :) < 0.2; %h2o conc are saturated at 30 mmol/mol 
bad_H2Oc_dl = bad_H2Oc_dl | D(78, :) < 0; %no variabiltiy 
D(86, bad_H2Oc_dl) = NaN;

%convert to mixing ratio
mixing_dl=D(86, :)./(1-(D(86, :)/1000)); %mol h2o/mol dry air
D(86, :)= mixing_dl; %mmol h2o/mol dry air

%clean GOES
H2Oc_goes_isnan = isnan(D(212,:))==1 | isempty(D(212,:)) ==1 | isfinite(D(212,:)) == 0;
bad_H2Oc_goes = H2Oc_goes_isnan ==1;
bad_H2Oc_goes = bad_H2Oc_goes | D(212, :) > 29.8; %h2o conc are saturated at 30 mmol/mol 
bad_H2Oc_goes = bad_H2Oc_goes | D(212, :) < 0.2; %h2o conc are saturated at 30 mmol/mol 
bad_H2Oc_goes = bad_H2Oc_goes | D(204, :) == 0; %no variabilty 
D(212, bad_H2Oc_goes) = NaN;

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
Tkact_mat   = (D(8,:)+273.15)  ./ (1 + 0.00032 * h2o_mol_frac); %[K]
Tkact_dl   = (D(84,:)+273.15)  ./ (1 + 0.00032 * h2o_mol_frac); %[K]
Tkact_goes   = (D(210,:)+273.15)  ./ (1 + 0.00032 * h2o_mol_frac); %[K]

Tkact_mat   = Tkact_mat - 273.15; %[C]
Tkact_dl   = Tkact_dl - 273.15; %[C]
Tkact_goes   = Tkact_goes - 273.15; %[C]
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
%37% of the observation will be rho_fixed at LowDes
%----------------------------------------------------------------------
%winds
%----------------------------------------------------------------------
%If wind vectors are wrong, then all wind and momentum fluxes are wrong
%clean mat
%component wind problems - variance is 0 is impossible
%u_wind
u_mat_var=D(9,:);
Umat_isnan=isnan(D(5,:))==1 | isempty(D(5,:)) ==1 | isfinite(D(5,:)) == 0;
ubad_wind_mat=u_mat_var(1,:)==0;
ubad_wind_mat = ubad_wind_mat | Umat_isnan(1,:)==1;
D(5,ubad_wind_mat)=NaN;

%v_wind
v_mat_var=D(10,:);
Vmat_isnan=isnan(D(6,:))==1 | isempty(D(6,:)) ==1 | isfinite(D(6,:)) == 0;
vbad_wind_mat = v_mat_var(1,:)==0;
vbad_wind_mat = vbad_wind_mat | Vmat_isnan(1,:)==1;
D(6,vbad_wind_mat)=NaN;

w_mat_var=D(11,:);
Wmat_isnan=isnan(D(7,:))==1 | isempty(D(7,:)) ==1 | isfinite(D(7,:)) == 0;
wbad_wind_mat = w_mat_var(1,:)==0;
wbad_wind_mat = wbad_wind_mat | Wmat_isnan(1,:)==1;
D(7,wbad_wind_mat)=NaN;


%clean dl
%DATA LOGGER
u_dl_var=D(66,:);%'Ux variance
Udl_isnan=isnan(D(81,:))==1 | isempty(D(81,:)) ==1 | isfinite(D(81,:)) == 0;
ubad_wind_dl=u_dl_var(1,:)==0;
ubad_wind_dl = ubad_wind_dl | Udl_isnan(1,:)==1;
D(81,ubad_wind_dl)=NaN;

v_dl_var=D(71,:);%'Uy variance
Vdl_isnan=isnan(D(82,:))==1 | isempty(D(82,:)) ==1 | isfinite(D(82,:)) == 0;
vbad_wind_dl = v_dl_var(1,:)==0;
vbad_wind_dl = vbad_wind_dl | Vdl_isnan(1,:)==1;
D(82,vbad_wind_dl)=NaN;

w_dl_var=D(60,:);% Uz variance
Wdl_isnan=isnan(D(83,:))==1 | isempty(D(83,:)) ==1 | isfinite(D(83,:)) == 0;
wbad_wind_dl = w_dl_var(1,:)==0;
wbad_wind_dl = wbad_wind_dl | Wdl_isnan(1,:)==1;
D(83,wbad_wind_dl)=NaN;


%clean goes
u_goes_var = D(192,:);% Ux variance
Ugoes_isnan=isnan(D(207,:))==1 | isempty(D(207,:)) ==1 | isfinite(D(207,:)) == 0;
ubad_wind_goes = u_goes_var(1,:) == 0;
ubad_wind_goes = ubad_wind_goes | Ugoes_isnan ==1;
ubad_wind_goes = ubad_wind_goes | D(207,:) < -50;
D(207,ubad_wind_goes) = NaN;% Ux

v_goes_var = D(197,:);% Uy variance 
Vgoes_isnan=isnan(D(208,:))==1 | isempty(D(208,:)) ==1 | isfinite(D(208,:)) == 0;
vbad_wind_goes = v_goes_var==0;
vbad_wind_goes = vbad_wind_goes | Vgoes_isnan ==1;
vbad_wind_goes = vbad_wind_goes | D(208,:) < -50;
D(208,vbad_wind_goes) = NaN;% Uy

w_goes_var = D(186,:);% Uz variance
Wgoes_isnan=isnan(D(209,:))==1 | isempty(D(209,:)) ==1 | isfinite(D(209,:)) == 0;
wbad_wind_goes = w_goes_var==0;
wbad_wind_goes = wbad_wind_goes | Wgoes_isnan ==1;
wbad_wind_goes = wbad_wind_goes | D(209,:) < -4000; % there are some -8000s in the dataset 
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
bad = bad | isnan(D(11,:)) == 1; %missing h2o
D(21,bad)=NaN;

%clean data logger
%bad = bad_wind_dl;
bad = D(60,:) == 0;
bad = bad | isnan(D(60,:)) == 1; % is there uz uz data
bad = bad | D(66,:) == 0; % sensor not working
bad = bad | isnan(D(66,:)) == 1; % is there ux ux data
D(61,bad)=NaN;

bad = D(60,:) == 0;
bad = bad | isnan(D(60,:)) == 1; % is there uz uz data
bad = bad | D(71,:) == 0; % sensor not working
bad = bad | isnan(D(71,:)) == 1; % is there uy uy data
D(62,bad)=NaN;

%clean GOES
bad = D(186,:) == 0; % uz uz
bad = bad | isnan(D(186,:)) == 1; % is there uz uz data
bad = bad | D(192,:) == 0; % sensor not working
bad = bad | isnan(D(192,:)) == 1; % is there ux ux data
bad = bad | D(187,:) <-700; % is there ux ux data
D(187,bad)=NaN;

bad = D(186,:) == 0; % uz uz
bad = bad | isnan(D(186,:)) == 1; % is there uz uz data
bad = bad | D(197,:) == 0; % sensor not working
bad = bad | isnan(D(197,:)) == 1; % is there uy uy data
bad = bad | D(197,:) <-700;
D(188,bad)=NaN;

%----------------------------------------------------------------------
%Sensible heat
%----------------------------------------------------------------------
%clean fast - this is despiked in fastflux
bad = D(12,:) == 0; % Ts var == 0 
bad = bad | isnan(D(12,:)) == 1; %missing Ts
bad = bad | D(11,:) == 0; %w variability = 0
bad = bad | isnan(D(11,:)) == 1; %missing w
bad = bad | time > 202.6 & time < 202.61;
D(29,bad)=NaN;

bad = D(29,:)<-350;
bad = bad | D(29,:)> 550;
D(29,bad)=NaN;

%clean data logger
%need a bad covariance term
bad = D(60,:) == 0;
bad = bad | isnan(D(60,:)) == 1; % is there uz uz data
bad = bad | isnan(D(80,:)) == 1; % is there uy uy data
bad = bad | D(80,:) == 0; % sensor not working
bad = bad | time > 202.6 & time < 202.61;
D(65,bad)=NaN;


%clean GOES 
bad = D(186,:) == 0; % uz uz
bad = bad | isnan(D(186,:)) == 1; % is there uz uz data
bad = bad | isnan(D(206,:)) == 1; % is there ux ux data
bad = bad | D(206,:) == 0; % sensor not working
bad = bad | time > 202.6 & time < 202.61;
bad = bad | D(191, :) < -50;
D(191, bad) = NaN;

%----------------------------------------------------------------------
%CO2 concentration
%----------------------------------------------------------------------
%clean fast
bad_fastCO2 = D(35,:) > 485;
bad_fastCO2 = bad_fastCO2 | D(35,:) < 380;
bad_fastCO2 = bad_fastCO2 | time > 1569.975 & time < 1569.98;
bad_fastCO2 = bad_fastCO2 | time > 1717.1 & time < 1717.33;
bad_fastCO2 = bad_fastCO2 | time > 1724.09 & time < 1724.24;
bad_fastCO2 = bad_fastCO2 | time > 1731.09 & time < 1731.4;
bad_fastCO2 = bad_fastCO2 | time > 382.97 & time < 383.01;
bad_fastCO2 = bad_fastCO2 | time > 332.745 & time < 332.8;
bad_fastCO2 = bad_fastCO2 | time > 678 & time < 680;
D(35,bad_fastCO2)=NaN;

%clean data logger
bad_dlCO2 = D(85,:) > 485;
bad_dlCO2 = bad_dlCO2 | D(85,:) < 380;
focus = time > 1710 & time < 1760;
bad_dlCO2 = bad_dlCO2 | focus & D(85,:) < 401;
bad_dlCO2 = bad_dlCO2 | time > 678 & time < 680;
D(85, bad_dlCO2) = NaN;

D(85,:)=D(85,:)./(1-mol_frac_dl); %co2/mol dry air

%clean GOES
bad_goesCO2 = D(211,:) > 485;
bad_goesCO2 = bad_goesCO2 | D(211,:) < 380;
focus = time > 1710 & time < 1760;
bad_goesCO2 = bad_goesCO2 | focus & D(211,:) < 401;
bad_goesCO2 = bad_goesCO2 | time > 678 & time < 680;
focus = time > 2150 & time < 2600;
bad_goesCO2 = bad_goesCO2 | focus & D(211,:) < 401;
D(211, bad_goesCO2) = NaN;

D(211,:)=D(211,:)./(1-mol_frac_goes); %co2/mol dry air
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
bad = bad | D(46,:) > 20;
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
bad_fco2_goes = bad_goesCO2;
bad_fco2_goes = bad_fco2_goes | isnan(D(186,:))==1; % sonic is not working no uz uz
bad_fco2_goes = bad_fco2_goes | D(186,:)==0;%no wind variability
bad_fco2_goes = bad_fco2_goes | D(201,:)==0; % co2 variability
bad_fco2_goes = bad_fco2_goes | isnan(D(201,:))==1; % sonic is not working no uz uz
bad_fco2_goes = bad_fco2_goes | isnan(D(189,:)) == 1; % no data
bad_fco2_goes = bad_fco2_goes | D(189,:) < -7000; 
D(189, bad_fco2_goes) = NaN;

%----------------------------------------------------------------------
%Adjust water fluxes using IRGA H2O variabilty compared with HMP variability
%NOTE: IRGA H2O variabilty seems high compared with HMP variability.
%Therefore, we over estimate FH2O.  This will reduce the water flux.

%Make sure the water vapor concentrations are cleaned first.
%----------------------------------------------------------------------
%IRGA calibrations
%The following is the data from the SiteDetails xls spreadsheet.  
Cnums=datenum(IRGAchanged(:,3), IRGAchanged(:,1), IRGAchanged(:,2));
Cnumsince2006=Cnums-datenum(2006,1,0);

firstday = min(time);
lastday = max(time);

length_Cnumsince2006=length(Cnumsince2006(:,1));
Cnumsince2006(1,1)=firstday-1;
Cnumsince2006(length_Cnumsince2006+1,1)=lastday;
length_Cnumsince2006 = length_Cnumsince2006+1;

% serviced figure
f1=(Cnumsince2006./Cnumsince2006)*40;
f2=(Cnumsince2006./Cnumsince2006)*-5;
f3=(Cnumsince2006./Cnumsince2006)*40;
F=[f1,f2,f3]';
[m,n]=size(F);
l=m*n;
F=reshape(F,l,1);
C06=[Cnumsince2006; Cnumsince2006; Cnumsince2006];
C06=sort(C06);

figure(2)
plot(time, D(100,:), 'k.');
hold on
plot(time, D(42,:), 'b.');
hold on
plot(time, D(126,:), 'g.');
hold on
plot(C06, F, 'r-');
legend('black=hmp', 'blue=irga', 'green=h2o zero (using span gas', 'irga calibration');

Pmat=[];
Pdl=[];
Pgoes=[];


for i=2:length_Cnumsince2006
    %proc data
    good_mat = time(1,:) < Cnumsince2006(i)-2 & time(1,:) > Cnumsince2006(i-1)+2 & isfinite(D(42,:))==1 & isfinite(D(100,:))==1;
    irga=D(42,good_mat);
    hmp=D(100,good_mat);
    p = polyfit(hmp,irga,1) %p(1) is slope; p(2) is intercept
    c_tm = time(1,:) <= Cnumsince2006(i) & time(1,:) > Cnumsince2006(i-1); %adjustment period
    
    %record so you can look later
    tt= median(time(1,good_mat)) 
    Pmat=[Pmat; tt, p];
    
    %The offset does not make sense here. Adjust water fluxes by dividing by the slope.
    outrageous_p = isnan(p(1,1))==1 | isfinite(p(1,1)) == 0 | p(1) > 1.3 | p(1) < 0.6
    
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
    good_dl = time(1,:) < Cnumsince2006(i)-2 & time(1,:) > Cnumsince2006(i-1)+2 & isfinite(D(86,:))==1 & isfinite(D(100,:))==1;
    irga=D(86,good_dl);
    hmp=D(100,good_dl);
    p = polyfit(hmp,irga,1) %p(1) is slope; p(2) is intercept
    c_tm = time(1,:) <= Cnumsince2006(i) & time(1,:) > Cnumsince2006(i-1); %adjustment period
    
    %record so you can look later
    tt= median(time(1,good_dl)) 
    Pdl=[Pdl; tt, p];
    
    outrageous_p = isnan(p(1,1))==1 | isfinite(p(1,1)) == 0 | p(1,1) > 1.3 | p(1,1) < 0.6
    
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
    good_goes = time(1,:) < Cnumsince2006(i)-2 & time(1,:) > Cnumsince2006(i-1)+2 & isfinite(D(212,:))==1 & isfinite(D(226,:))==1;
    irga=D(212,good_goes);
    hmp=D(226,good_goes);
    p = polyfit(hmp,irga,1) %p(1) is slope; p(2) is intercept
    c_tm = time(1,:) <= Cnumsince2006(i) & time(1,:) > Cnumsince2006(i-1); %adjustment period
    
    %record so you can look later
    tt= median(time(1,good_goes)) 
    Pgoes=[Pgoes; tt, p];
    
    %The offset does not make sense here. Adjust water fluxes by dividing by the slope.
    outrageous_p = isnan(p(1,1))==1 | isfinite(p(1,1)) == 0 | p(1,1) > 1.3 | p(1,1) < 0.6
    
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
save('C:\towerData\combined\CombinedInfo\DC_LowDes\HMP_h2o_calibration_DC_LowDes.mat', 'Pmat', 'Pdl', 'Pgoes');
%----------------------------------------------------------------------
%FH2O
%----------------------------------------------------------------------
%clean fast

bad_mat_et = bad_H2Oc_mat;
bad_mat_et = bad_mat_et | D(43,:) == 0; %
bad_mat_et = bad_mat_et | isnan(D(43,:)) == 1; %
bad_mat_et = bad_mat_et | D(11,:) == 0; %uz variability = 0
bad_mat_et = bad_mat_et | isnan(D(11,:)) == 1; %missing uz
bad_mat_et = bad_mat_et | D(47,:) < -8;
bad_mat_et = bad_mat_et | D(47,:) > 20;

D(47, bad_mat_et) = NaN;

%clean data logger
bad_dl_et = bad_H2Oc_dl;
bad_dl_et = bad_dl_et | isnan(D(60,:))==1; % sonic is not working no uz uz
bad_dl_et = bad_dl_et | D(60,:)==0;%no wind variability
bad_dl_et = bad_dl_et | D(78,:)==0; % h2o variability
bad_dl_et = bad_dl_et | isnan(D(78,:))==1; % sonic is not working no uz uz
bad_dl_et = bad_dl_et | isnan(D(64,:)) == 1; % no data

D(64, bad_dl_et) = NaN;

%clean GOES
bad_goes_et = bad_H2Oc_goes;
bad_goes_et = bad_goes_et | isnan(D(186,:))==1; % sonic is not working no uz uz
bad_goes_et = bad_goes_et | D(186,:)==0;%no wind variability
bad_goes_et = bad_goes_et | D(204,:)==0; % h2o variability
bad_goes_et = bad_goes_et | isnan(D(204,:))==1; % 
bad_goes_et = bad_goes_et | isnan(D(190,:)) == 1; % no data
bad_goes_et = bad_goes_et | D(190,:) < -6000;
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
bad = D(87,:) > 760; 
bad = bad | D(87,:) < -400;
bad = bad | time > 170 & time < 332.8; % check this out...
bad = bad | time > 451.84 & time < 452.785;
bad = bad | time > 521.5 & time < 552.95;
bad = bad | time > 611.94 & time < 613.12;
bad = bad | time > 1297.02 & time < 1297.05;
D(87,bad) = NaN;

%clean GOES
bad = D(213,:) > 760; 
bad = bad | D(213,:) < -400;
bad = bad | time > 170 & time < 332.8; % check this out...
bad = bad | time > 451.84 & time < 452.785;
bad = bad | time > 521.5 & time < 552.95;
bad = bad | time > 611.94 & time < 613.12;
bad = bad | time > 1297.02 & time < 1297.05;
D(213,bad) = NaN;

%----------------------------------------------------------------------
%Incoming solar radiation - pyranometer = SOLAR_IN
%----------------------------------------------------------------------
%looks like pyrr in and out are flipped ==> calibrations are wrong
In_dl=D(89,:);
Out_dl=D(88,:);

In_goes=D(215,:);
Out_goes=D(214,:);
%adjust the calibration
%ex. curr_multiplier = 1e-3*1/CFACTOR1(good); 
%D(89,:)= D(89,:).* curr_multiplier; %dl
curr_m_dn = 1e-3*1/0.00002246; % see SiteDetails.xlsx
curr_m_up=1e-3*1/0.00001969;
  
D(88,:)=In_dl(1,:)*(curr_m_up/curr_m_dn);
 
bad = D(88,:) < -10;
bad = bad | time > 824.86 & time < 837.2; 
D(88,bad) = NaN;

%clean goes
D(214,:)= In_goes(1,:)*(curr_m_up/curr_m_dn);

bad = D(214,:) < -10;
bad = bad | time > 824.86 & time < 837.2; 
bad = bad | time > 452.72 & time < 452.78;
bad = bad |  time < 119; 
D(214,bad) = NaN;

%----------------------------------------------------------------------
%Outgoing solar radiation - pyranometer = SOLAR_OUT
%----------------------------------------------------------------------
%clean data logger
D(89,:) = Out_dl;
bad = D(89,:) < -100;
D(89,bad)= NaN;

%clean goes
D(215,:) = Out_goes;
bad = D(215,:) < -100;
D(215,bad)= NaN;
%----------------------------------------------------------------------
%PAR_In
%----------------------------------------------------------------------
%clean data logger
bad = D(90,:) < -20; 
D(90,bad) = NaN;

%clean GOES
bad = D(216,:) < -20;
D(216,bad) = NaN;

%----------------------------------------------------------------------
%PAR_Out
%----------------------------------------------------------------------
%clean data logger
bad = D(91,:) <-10; 
bad = bad | time > 655.85 & time < 700.5; 
bad = bad | time > 1136.5 & time < 1159.872; 
bad = bad | time > 1177.76 & time < 1287.995; 
bad = bad | time > 1746.74 & time < 1749.956; 
bad = bad | time > 1797.84 & time < 1837.95; 
bad = bad | time > 1923.688 & time < 1936; 
D(91,bad) = NaN;

%clean GOES
bad = D(217,:) <-10; 
bad = bad | time > 655.85 & time < 700.5; 
bad = bad | time > 836.82 & time < 836.905; 
bad = bad | time > 1136.5 & time < 1159.872; 
bad = bad | time > 1177.76 & time < 1287.995; 
bad = bad | time > 1746.74 & time < 1749.956; 
bad = bad | time > 1797.84 & time < 1837.95; 
bad = bad | time > 1923.688 & time < 1936; 
bad = bad | time > 789.745 & time < 789.775; 
bad = bad | time > 1580.76 & time < 1580.78; 
bad = bad | time > 1171.77 & time < 1171.78;
D(217,bad) = NaN;
%----------------------------------------------------------------------
%T_107
%----------------------------------------------------------------------
%dl
bad=D(102,:)<-20;
D(102,bad) = NaN;

%goes
bad=D(228,:)<0.2;
D(228,bad) = NaN;

%----------------------------------------------------------------------
%Rain
%----------------------------------------------------------------------
%data logger
bad = D(103,:)>20;
D(103,bad) = NaN;

%goes
bad = D(229,:)>20;
D(229,bad) = NaN;

%----------------------------------------------------------------------
%NDVI
%----------------------------------------------------------------------
bad = D(96,:) > 1;
D(96,bad)=NaN;

bad = D(222,:) > 1;
D(222,bad)=NaN;

%----------------------------------------------------------------------
%Soil Moisture
%----------------------------------------------------------------------
%remove soil moisture sensor 1
D(282,:) = NaN;

%----------------------------------------------------------------------
%Soil T
%----------------------------------------------------------------------


%----------------------------------------------------------------------
%LWS
%----------------------------------------------------------------------


%----------------------------------------------------------------------
%Matric
%----------------------------------------------------------------------
%no sensors

elseif a==2
  
%Mop-up

TEMP =[];
H2O =[];
rho_fix =[];
    
  
%Sensible Heat
%dl
bad = D(339,:)<-350;
bad = bad | D(339,:)> 550;
D(339,bad)=NaN;

%goes
bad = D(340,:)<-350;
bad = bad | D(340,:)> 550;
D(340, bad) = NaN;

%Fh2o
%dl
bad_dl_et = D(347,:) < -8;
bad_dl_et = bad_dl_et | D(347,:) > 12;
D(347, bad_dl_et) = NaN;

%goes
bad_goes_et = D(348,:) < -8;
bad_goes_et = bad_goes_et |  D(348,:) > 12;
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
bad = bad | D(345,:) > 20;
D(345, bad) = NaN;

%goes
bad = D(346,:) < -20;
bad = bad | D(346,:) > 20;
D(346, bad) = NaN;

%Ustar
%clean GOES 
bad = D(323,:) > 100; 
D(323,bad)=NaN;

end