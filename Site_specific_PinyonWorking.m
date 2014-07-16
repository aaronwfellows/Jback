function [D, TEMP, H2O, rho_fix, v_wind] = Site_specific_PinyonWorking(HEADER, D, Press, R_mol, Tc, a)

%This file is a site specific clean-up for DMERGE data at the DC_Pinyon
%awf 10/23/2012

%OUTPUTS:
%D is the filtered DMERGE data file
%TEMP is a cleaned and filled temperature in deg C 
%H2O is the water mixing ratio [mmol h2o/ mol of dry air]
%rho_fix is bad temp and water vapor measurements
%----> rho_fix is a handle that points to fast calculated
%----> fluxes that used bad dry density estimates 
%v_wind is cleaned up bad wind components

%When was the IRGA serviced?
IRGAchanged = [5,17,2006;1,10,2007;7,5,2007;2,27,2008;8,21,2008;2,26,2009;3,19,2010;3,16,2011;8,9,2011;6,12,2012;10,2,2013];

if a==1
%----------------------------------------------------------------------
%fix up the goes wind problem in the u direction
focus= D(1,:) > 529.26 & D(1,:) < 535.8955;
focus= focus | D(1,:) > 1266.065;
D(207,~focus)= D(207,~focus).*100; %this is divided  
%----------------------------------------------------------------------

%use ones through the program
len=length(D(1,:));
o=ones(1,len);

%use time to remove bad data
time = D(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%By measurement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------
%HMP_Temp 
%----------------------------------------------------------------------
%data logger
badHMP_dl = isnan(D(97,:))==1 | isempty(D(97,:)) ==1 | isfinite(D(97,:)) == 0;
badHMP_dl = badHMP_dl | D(97,:) < -15;
badHMP_dl = badHMP_dl | time > 501.57 & time < 508.1;
badHMP_dl = badHMP_dl | time > 586.08 & time < 586.12;
badHMP_dl = badHMP_dl | time > 645.442 & time < 646.4;
badHMP_dl = badHMP_dl | time > 649.76 & time < 650.22;
badHMP_dl = badHMP_dl | time > 650.8 & time < 651.15;
badHMP_dl = badHMP_dl | time > 652.78 & time < 653.32;
badHMP_dl = badHMP_dl | time > 653.76 & time < 654.56;
badHMP_dl = badHMP_dl | time > 654.86 & time < 655.02;
badHMP_dl = badHMP_dl | time > 657.946 & time < 658.2;
badHMP_dl = badHMP_dl | time > 661.81 & time < 679.57;
badHMP_dl = badHMP_dl | time > 710 & time < 711.95;
badHMP_dl = badHMP_dl | time > 1536.1 & time < 1539.89;
D(97,badHMP_dl) = NaN;

%goes
badHMP_goes = isnan(D(223,:))==1 | isempty(D(223,:)) ==1 | isfinite(D(223,:)) == 0;
badHMP_goes = badHMP_goes | D(223,:) < -8;
badHMP_goes = badHMP_goes  | time > 501.57 & time < 508.1;
badHMP_goes = badHMP_goes | time > 586.08 & time < 586.12;
badHMP_goes = badHMP_goes | time > 645.442 & time < 646.4;
badHMP_goes = badHMP_goes | time > 649.76 & time < 650.22;
badHMP_goes = badHMP_goes | time > 650.8 & time < 651.15;
badHMP_goes = badHMP_goes | time > 652.78 & time < 653.32;
badHMP_goes = badHMP_goes | time > 653.76 & time < 654.56;
badHMP_goes = badHMP_goes | time > 654.86 & time < 655.02;
badHMP_goes = badHMP_goes | time > 657.946 & time < 658.2;
badHMP_goes = badHMP_goes | time > 661.81 & time < 679.57;
badHMP_goes = badHMP_goes  | time > 1536.1 & time < 1539.89;
badHMP_goes = badHMP_goes | time > 1110.88 & time < 1446.791; % time shifted
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
bad_Ts_mat= Ts_mat_isnan == 1 | time>1574.0208 & time <1574.0209;
bad_Ts_mat= bad_Ts_mat | time>654.999 & time <655.001;
bad_Ts_mat= bad_Ts_mat | time>174.8958 & time <174.8959;
D(8,bad_Ts_mat)=NaN;

%clean dl
Ts_dl_isnan=isnan(D(84,:))==1 | isempty(D(84,:)) ==1 | isfinite(D(84,:)) == 0;
bad_Ts_dl= Ts_dl_isnan ==1;
D(84,bad_Ts_dl)=NaN;

%clean GOES
Ts_goes_isnan=isnan(D(210,:))==1 | isempty(D(210,:)) ==1 | isfinite(D(210,:)) == 0;
bad_Ts_goes = Ts_goes_isnan ==1 | D(210,:)==0;
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
bad = time > 436.74 & time < 446.965;
bad = bad | time > 388.82 & time < 389.22;
bad = bad | time > 395.25 & time < 396.99;
bad = bad | time > 429.8 & time < 430.3;
bad = bad | time > 586.1 & time < 655.03;
bad = bad | time > 704 & time < 707;
bad = bad | time > 710.7 & time < 712;
bad = bad | time > 734 & time < 737.5;
bad = bad | time > 767.86 & time < 789.1;
bad = bad | time > 895 & time < 896.72;
bad = bad | time > 943.24 & time < 964.8;
bad = bad | time > 1085 & time < 1154.675;
bad = bad | time > 1901.2 & time < 1987.885;
bad = bad | time > 2861.5 & time < 2877.2;
focus = time > 0 & time < 934 | time > 1000;
bad = bad | focus & D(98,:) < 0.01;
D(98,bad) = NaN;

%goes
bad = time > 436.74 & time < 446.965;
bad = bad | time > 388.82 & time < 389.22;
bad = bad | time > 395.25 & time < 396.99;
bad = bad | time > 429.8 & time < 430.3;
bad = bad | time > 586.1 & time < 655.03;
bad = bad | time > 704 & time < 707;
bad = bad | time > 710.7 & time < 712;
bad = bad | time > 734 & time < 737.5;
bad = bad | time > 767.86 & time < 789.1;
bad = bad | time > 895 & time < 896.72;
bad = bad | time > 943.24 & time < 964.8;
bad = bad | time > 1085 & time < 1154.675;
bad = bad | time > 1266.8 & time < 1269.5;
bad = bad | time > 1901.2 & time < 1987.885;
bad = bad | time > 1110.88 & time < 1446.791; % time shifted - could get some of this data back
bad = bad | time > 2861.5 & time < 2877.2;
bad = bad | time > 2992 & time < 2999;
focus = (time > 0 & time < 934) | (time > 1000 & time < 2650)| (time > 2880 & time < 2920);
bad = bad | focus & D(224,:) < 0.01;
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
bad_H2Oc_mat = bad_H2Oc_mat |D(42,:) < 0;
bad_H2Oc_mat = bad_H2Oc_mat |D(42,:) > 29.9;
bad_H2Oc_mat = bad_H2Oc_mat | time > 1088.6 & time < 1138.5; % compared with DC_Burn
bad_H2Oc_mat = bad_H2Oc_mat | time > 1453 & time < 1462;
bad_H2Oc_mat = bad_H2Oc_mat | time > 1503 & time < 1529;
D(42, bad_H2Oc_mat) = NaN;

%clean data logger
H2Oc_dl_isnan = isnan(D(86,:))==1 | isempty(D(86,:)) ==1 | isfinite(D(86,:)) == 0;
bad_H2Oc_dl = H2Oc_dl_isnan ==1;
bad_H2Oc_dl = bad_H2Oc_dl | D(86,:) < 0;
bad_H2Oc_dl = bad_H2Oc_dl | D(86,:) > 29.9;
bad_H2Oc_dl = bad_H2Oc_dl | time > 1088.6 & time < 1138.5;
bad_H2Oc_dl = bad_H2Oc_dl | time > 1453 & time < 1462;
bad_H2Oc_dl = bad_H2Oc_dl | time > 1503 & time < 1529;
bad_H2Oc_dl = bad_H2Oc_dl | time > 2981.72 & time < 2992.04;
D(86, bad_H2Oc_dl) = NaN;

%convert to mixing ratio
mixing_dl=D(86, :)./(1-(D(86, :)/1000)); %mol h2o/mol dry air
D(86, :)= mixing_dl; %mmol h2o/mol dry air

%clean GOES
%clean GOES
H2Oc_goes_isnan = isnan(D(212,:))==1 | isempty(D(212,:)) ==1 | isfinite(D(212,:)) == 0;
bad_H2Oc_goes = H2Oc_goes_isnan ==1;
bad_H2Oc_goes = bad_H2Oc_goes | D(212,:) < 0;
bad_H2Oc_goes = bad_H2Oc_goes | D(212,:) > 29.9;
bad_H2Oc_goes  = bad_H2Oc_goes  | time > 1088.6 & time < 1138.5;
bad_H2Oc_goes = bad_H2Oc_goes | time > 1453 & time < 1462;
bad_H2Oc_goes = bad_H2Oc_goes | time > 1503 & time < 1529;
bad_H2Oc_goes = bad_H2Oc_goes | time > 2981.72 & time < 2992.04;
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
%~16% will be rho_fixes
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
D(207,ubad_wind_goes) = NaN;% Ux

v_goes_var = D(197,:);% Uy variance 
Vgoes_isnan=isnan(D(208,:))==1 | isempty(D(208,:)) ==1 | isfinite(D(208,:)) == 0;
vbad_wind_goes = v_goes_var==0;
vbad_wind_goes = vbad_wind_goes | Vgoes_isnan ==1;
D(208,vbad_wind_goes) = NaN;% Uy

w_goes_var = D(186,:);% Uz variance
Wgoes_isnan=isnan(D(209,:))==1 | isempty(D(209,:)) ==1 | isfinite(D(209,:)) == 0;
wbad_wind_goes = w_goes_var==0;
wbad_wind_goes = wbad_wind_goes | Wgoes_isnan ==1;
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
D(187,bad)=NaN;

bad = D(186,:) == 0; % uz uz
bad = bad | isnan(D(186,:)) == 1; % is there uz uz data
bad = bad | D(197,:) == 0; % sensor not working
bad = bad | isnan(D(197,:)) == 1; % is there uy uy data
D(188,bad)=NaN;

%----------------------------------------------------------------------
%Sensible heat
%----------------------------------------------------------------------
%clean fast - this is despiked in fastflux
bad = D(12,:) == 0; % Ts var == 0 
bad = bad | isnan(D(12,:)) == 1; %missing Ts
bad = bad | D(11,:) == 0; %w variability = 0
bad = bad | isnan(D(11,:)) == 1; %missing w
D(29,bad)=NaN;

bad = D(29,:)<-350;
bad = bad | D(29,:)> 700;
D(29,bad)=NaN;

%clean data logger
%need a bad covariance term
bad = D(60,:) == 0;
bad = bad | isnan(D(60,:)) == 1; % is there uz uz data
bad = bad | isnan(D(80,:)) == 1; % is there uy uy data
bad = bad | D(80,:) == 0; % sensor not working
D(65,bad)=NaN;


%clean GOES 
bad = D(186,:) == 0; % uz uz
bad = bad | isnan(D(186,:)) == 1; % is there uz uz data
bad = bad | isnan(D(206,:)) == 1; % is there ux ux data
bad = bad | D(206,:) == 0; % sensor not working
D(191, bad) = NaN;

%----------------------------------------------------------------------
%CO2 concentration
%----------------------------------------------------------------------
%clean fast
bad_fastCO2 = D(35,:) > 476;
bad_fastCO2 = bad_fastCO2 | D(35,:) < 360;
focus = time > 2660 & time < 2700;
bad_fastCO2 = bad_fastCO2 | (focus & D(85,:) < 390);
D(35,bad_fastCO2)=NaN;

%clean data logger
bad_dlCO2 = D(85,:) > 470;
bad_dlCO2 = bad_dlCO2 | D(85,:) < 360;
focus = time > 1902.5;
bad_dlCO2 = bad_dlCO2 | focus & D(85,:) < 380.8;
focus = time > 2100;
bad_dlCO2 = bad_dlCO2 | focus & D(85,:) < 388;
D(85, bad_dlCO2) = NaN;

D(85,:)=D(85,:)./(1-mol_frac_dl); %co2/mol dry air

%clean GOES
bad_goesCO2 = D(211,:) > 470;
bad_goesCO2 = bad_goesCO2 | D(211,:) < 360;
focus = time > 1902.5;
bad_goesCO2 = bad_goesCO2| focus & D(211,:) < 380.8;
focus = time > 2100;
bad_goesCO2 = bad_goesCO2 | focus & D(211,:) < 388;
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
D(189, bad_fco2_goes) = NaN;

%----------------------------------------------------------------------
%Adjust water fluxes using IRGA H2O variabilty compared with HMP variability
%NOTE: IRGA H2O variabilty seems high compared with HMP variability.
%Therefore, we over estimate FH2O.  This will reduce the water flux.

%Make sure the water vapor concentrations are cleaned first.
%----------------------------------------------------------------------
%IRGA calibrations
%The following is the data from the SiteDetails xls spreadsheet - from Greg's notebook.  
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
legend('black=hmp', 'blue=irga', 'irga calibration');

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
    outrageous_p = isnan(p(1,1))==1 | isfinite(p(1,1)) == 0 | p(1,1) > 1.5 | p(1,1) < 0.8
    
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
    
    outrageous_p = isnan(p(1,1))==1 | isfinite(p(1,1)) == 0 | p(1,1) > 1.5 | p(1,1) < 0.8
    
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
    outrageous_p = isnan(p(1,1))==1 | isfinite(p(1,1)) == 0 | p(1,1) > 1.5 | p(1,1) < 0.8
    
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

save('C:\towerData\combined\CombinedInfo\DC_Pinyon\HMP_h2o_calibration_DC_Pinyon.mat', 'Pmat', 'Pdl', 'Pgoes');
close(figure(2))
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
bad_mat_et = bad_mat_et | D(47,:) > 12;

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
bad = bad | time > 209.2 & time < 376;
bad = bad | time > 240 & time < 374.1;
bad = bad | time > 1453.5 & time < 1531.92; 
bad = bad | time > 1978.15 & time < 1985.1; 
bad = bad | time > 1990.1 & time < 2010.28;
D(87,bad) = NaN;

%clean GOES
bad = D(213,:) > 760; 
bad = bad | D(213,:) < -200; 
bad = bad | time > 209.2 & time < 376; 
bad = bad | time > 428.2 & time < 433.75; 
bad = bad | time > 434.72 & time < 434.76;
bad = bad | time > 434.98 & time < 435.06;
bad = bad | time > 1453.5 & time < 1531.92; 
bad = bad | time > 1978.15 & time < 1985.1; 
bad = bad | time > 1990.1 & time < 2010.28;
bad = bad | time > 2014.8 & time < 2048.03;
focus = time > 185 & time < 195;
bad = bad | focus & D(213,:)==0;
focus = time > 430 & time < 460;
bad = bad | focus & D(213,:)==0;
focus = time > 1539 & time < 1541;
bad = bad | focus & D(213,:)==0;
D(213,bad) = NaN;

%----------------------------------------------------------------------
%Incoming solar radiation - pyranometer = SOLAR_IN
%----------------------------------------------------------------------
%clean data logger

%clean goes
focus = time > 1560.75 & time < 1560.9;
bad = bad | focus & D(214,:)==0;
D(214,bad) = NaN;
%----------------------------------------------------------------------
%Outgoing solar radiation - pyranometer = SOLAR_OUT
%----------------------------------------------------------------------
%clean data logger
bad = D(89,:) < -10;
bad = bad | time > 1479.945 & time < 1531.905;
bad = bad | time > 1539.87 & time < 1560.89;
D(89,bad) = NaN;

%clean goes 
bad = D(215,:) < -10;
bad = bad | time > 1479.945 & time < 1531.905;
bad = bad | time > 1539.87 & time < 1560.89; 
D(215,bad) = NaN;

%----------------------------------------------------------------------
%PAR_In
%----------------------------------------------------------------------
%clean data logger
bad = D(90,:) < -20; 
D(90,bad) = NaN;

%clean GOES
bad = D(216,:) < -20;
focus = time > 1560.75 & time < 1560.9;
bad = bad | focus & D(216,:)==0;
D(216,bad) = NaN;

%----------------------------------------------------------------------
%PAR_Out
%----------------------------------------------------------------------
%clean data logger
bad = D(91,:) <-20; 
bad = bad | time > 2841.3 & time < 2842.3;
D(91,bad) = NaN;

%clean GOES
bad = time > 789.03 & time < 808.9; 
bad = bad | time > 2841.3 & time < 2842.3;
bad = bad | time > 2848.05 & time < 2876.93;
D(217,bad) = NaN;

%----------------------------------------------------------------------
%T_107
%----------------------------------------------------------------------
%dl
bad = D(102,:) < -20;
bad = bad | time > 1536.01 & time < 1539.89;
D(102,bad) = NaN;

%goes
bad = D(228,:) < -20;
bad = bad | time > 1536.01 & time < 1539.89;
focus = time > 200 & time < 1400;
bad = bad | focus & D(228,:)==0;
bad = bad | time > 1110.88 & time < 1446.791; % time shifted
D(228,bad) = NaN;

%----------------------------------------------------------------------
%Rain
%----------------------------------------------------------------------
%data logger
bad = D(103,:)>100;
D(103,bad) = NaN;

%goes
bad = D(229,:)>100;

bad = bad | time > 1110.88 & time < 1446.791; % time shifted
D(229,bad) = NaN;

%----------------------------------------------------------------------
%NDVI
%----------------------------------------------------------------------
bad = D(96,:) > 1;
D(96,bad)=NaN;

bad = D(222,:) > 1;
bad = bad | time > 1110.88 & time < 1446.791; % time shifted
D(222,bad)=NaN;
D(218,bad)=NaN; % Red  In
D(220,bad)=NaN; % Red  Out
D(219,bad)=NaN; % NIR
D(221,bad)=NaN; % NIR
%----------------------------------------------------------------------
%Soil T
%----------------------------------------------------------------------
%T1
bad = D(288,:) < -10;
bad = bad | time > 1256.76 & time < 1440; 
D(288,bad)=NaN;

%T2
bad = D(289,:) < -10;
D(289,bad)=NaN;

%T3
bad = D(290,:) < -10;
D(290,bad)=NaN;

%T4
bad = D(291,:) < -10;
D(291,bad)=NaN;

%----------------------------------------------------------------------
%Clean Soil Moisture
%----------------------------------------------------------------------

%----------------------------------------------------------------------
%Fuel_M
%----------------------------------------------------------------------
%data logger - these can be cleaned up but they are in bad shape 
%right now I will kill all the data 10/23/2012
%D(288,:)=NaN;
%D(289,:)=NaN;

%----------------------------------------------------------------------
%LWS
%----------------------------------------------------------------------
%don't really know what this data should look like
D(293,:)=NaN;
D(294,:)=NaN;
D(295,:)=NaN;

%----------------------------------------------------------------------
%Clean Matric potential sensors 
%----------------------------------------------------------------------
bad = time > 2832 & time < 2834;
D(298,bad)=NaN;

bad = time > 1107.9 & time < 1160.992;
bad = bad | time > 1004.86 & time < 1004.92;
D(301,bad)=NaN;

focus = time > 1107.9 & time < 1164;
bad = focus & D(307,:) < 13;
bad = bad | D(307,:) > 24;
D(307,bad)=NaN;

bad = time > 2489.21 & time < 2993;
D(314,bad) = NaN;

bad = time > 1004.86 & time < 1004.9;
bad = bad | time > 1161.08 & time < 1161.085;
bad = bad | time > 2832 & time < 2834;
bad = bad | time > 2982.26 & time < 2982.3;
D(316,bad)=NaN;

bad = time > 1004.86 & time < 1004.9;
bad = bad | time > 1153.63 & time < 1161.12;
bad = bad | time > 2982.26 & time < 2982.3;
D(317,bad)=NaN;


bad = time > 1004.86 & time < 1004.9;
bad = bad | time > 1107.95 & time < 1161.12;
bad = bad | time > 2982.26 & time < 2982.3;
D(318,bad)=NaN;
bad = bad | time > 1560.8 & time < 2993;
D(319,bad)=NaN;


elseif a==2
  
%Mop-up

TEMP =[];
H2O =[];
rho_fix =[];
v_wind =[];    
  
%Sensible Heat
%dl
bad = D(339,:)<-350;
bad = bad | D(339,:)> 700;
D(339,bad)=NaN;

%goes
bad = D(340,:)<-350;
bad = bad | D(340,:)> 700;
D(340, bad) = NaN;


%Fh2o
%dl
bad_dl_et = D(347,:) < -8;
bad_dl_et = bad_dl_et | D(347,:) > 12;
D(347, bad_dl_et) = NaN;

%goes
bad_goes_et = D(348,:) < -8;
bad_goes_et = bad_goes_et | D(348,:) > 12;
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
%clean dl

%clean GOES


end