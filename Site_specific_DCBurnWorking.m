function [D, TEMP, H2O, rho_fix, v_wind] = Site_specific_DCBurnWorking(HEADER, D, Press, R_mol, Tc, a)

%When were the calibration changes for the irga?
IRGAchanged =[5,17,2006;1,10,2007;5,7,2007;10,23,2007;2,27,2008;3,19,2008;2,26,2009;3,17,2010;4,29,2010;3,16,2011;8,9,2011;6,12,2012];

%sonic temp problem --> proc fast shows low variability ~ day 700 to day
%965 & day 126 to 132 compared with HMP==> this could influence the
%sensible heat flux...must track this down awf 10/26/2012

%CO2 concentration ~days 650 - 750 is low - particulatly evident in data
%logger and goes data streams (perhaps filtered out in fastflux)...I don't
%know why yet but should understand this one. There is also increased CO2 conc variabilty 
%and the pattern is evident in the FCO2...I think this is probably a measurement problem 
%but I don't understand what is happening yet. - awf 10/26/2012


%This file is a site specific clean-up for DMERGE data at the DC_Burn
%awf 10/23/2012

% co2 concentrations are messed up around 
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
focus= D(1,:) > 529.26 & D(1,:) < 535.88;
focus= focus | D(1,:) > 1265;
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
badHMP_dl = badHMP_dl | time > 132.68 & time < 133.012;
badHMP_dl = badHMP_dl | time > 174.6 & time < 174.821;
badHMP_dl = badHMP_dl | time > 207.8 & time < 208.01;
badHMP_dl = badHMP_dl | time > 220.91 & time < 220.97;
badHMP_dl = badHMP_dl | time > 282.88 & time < 282.96;
badHMP_dl = badHMP_dl | time > 332.7 & time < 333.005;
badHMP_dl = badHMP_dl | time > 382.975 & time < 383.005;
badHMP_dl = badHMP_dl | time > 452.7 & time < 452.86;
badHMP_dl = badHMP_dl | time > 515.8 & time < 516.01;
badHMP_dl = badHMP_dl | time > 554.306 & time < 555.002;
badHMP_dl = badHMP_dl | time > 585.94 & time < 585.98;
badHMP_dl = badHMP_dl | time > 655.72 & time < 656.01;
badHMP_dl = badHMP_dl | time > 789.71 & time < 790.01;
badHMP_dl = badHMP_dl | time > 1171.76 & time < 1172.01;
badHMP_dl = badHMP_dl | time > 1288.11 & time < 1289.01;
badHMP_dl = badHMP_dl | time > 1414.035 & time < 1415.01;
badHMP_dl = badHMP_dl | time > 1550.88 & time < 1551.01;
badHMP_dl = badHMP_dl | time > 1580.8958 & time < 1580.896;
badHMP_dl = badHMP_dl | time > 1585.88 & time < 1586.01;
badHMP_dl = badHMP_dl | time > 1746.7 & time < 1747.01;
badHMP_dl = badHMP_dl | time > 1749.905 & time < 1750.01;
badHMP_dl = badHMP_dl | time > 1838.05 & time < 1839.004;
badHMP_dl = badHMP_dl | time > 1861.055 & time < 1862.01;
badHMP_dl = badHMP_dl | time > 1887.95 & time < 1888.01;
badHMP_dl = badHMP_dl | D(97,:) < -35;
D(97,badHMP_dl) = NaN;

%goes
badHMP_goes = isnan(D(223,:))==1 | isempty(D(223,:)) ==1 | isfinite(D(223,:)) == 0;
badHMP_goes = badHMP_goes | D(223,:) < -30;
badHMP_goes = badHMP_goes | time > 585.94 & time < 585.98;
badHMP_goes = badHMP_goes | time > 1181.08 & time < 1181.09;
badHMP_goes = badHMP_goes | time > 1193.872 & time < 1193.876;
badHMP_goes = badHMP_goes | time > 1254.68 & time < 1254.7;
badHMP_goes = badHMP_goes | time > 1256.24 & time < 1256.26;
badHMP_goes = badHMP_goes | time > 1294.68 & time < 1294.7;
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
bad_Ts_mat= Ts_mat_isnan == 1 | time > 132.68 & time < 133.012;
bad_Ts_mat = bad_Ts_mat | time > 167.9582 & time < 167.9583;
bad_Ts_mat = bad_Ts_mat | time > 174.675 & time < 174.85;
bad_Ts_mat = bad_Ts_mat | time > 207.102 & time < 207.106;
bad_Ts_mat = bad_Ts_mat | time > 207.84 & time < 208.01;
bad_Ts_mat = bad_Ts_mat | time > 220.9 & time < 220.99;
bad_Ts_mat = bad_Ts_mat | time > 282.97 & time < 283.01;
bad_Ts_mat = bad_Ts_mat | time > 332.76 & time < 333.01;
bad_Ts_mat = bad_Ts_mat | time > 374.07 & time < 374.09;
bad_Ts_mat = bad_Ts_mat | time > 382.995 & time < 383.005;
bad_Ts_mat = bad_Ts_mat | time > 452.72 & time < 452.88;
bad_Ts_mat = bad_Ts_mat | time > 492.915 & time < 492.92;
bad_Ts_mat = bad_Ts_mat | time > 515.86 & time < 516.01;
bad_Ts_mat = bad_Ts_mat | time > 554.36 & time < 555.015;
bad_Ts_mat = bad_Ts_mat | time > 654.935 & time < 654.945;
bad_Ts_mat = bad_Ts_mat | time > 655.715 & time < 656.01;
bad_Ts_mat = bad_Ts_mat | time > 700.2 & time < 701.01;
bad_Ts_mat = bad_Ts_mat | time > 789.76 & time < 790.01;
D(8,bad_Ts_mat)=NaN;

%clean dl
Ts_dl_isnan=isnan(D(84,:))==1 | isempty(D(84,:)) ==1 | isfinite(D(84,:)) == 0;
bad_Ts_dl= Ts_dl_isnan ==1;
bad_Ts_dl = bad_Ts_dl | time > 132.68 & time < 133.012;
bad_Ts_dl = bad_Ts_dl | time > 174.6 & time < 174.85;
bad_Ts_dl = bad_Ts_dl | time > 207.102 & time < 207.106;
bad_Ts_dl = bad_Ts_dl | time > 207.84 & time < 208.01;
bad_Ts_dl = bad_Ts_dl | time > 220.9 & time < 220.99;
bad_Ts_dl = bad_Ts_dl | time > 282.89 & time < 282.97;
bad_Ts_dl = bad_Ts_dl | time > 332.69 & time < 333.015;
bad_Ts_dl = bad_Ts_dl | time > 382.96 & time < 383.01;
bad_Ts_dl = bad_Ts_dl | time > 452.7 & time < 452.88;
bad_Ts_dl = bad_Ts_dl | time > 515.8 & time < 516.01;
bad_Ts_dl = bad_Ts_dl | time > 554.3 & time < 555.015;
bad_Ts_dl = bad_Ts_dl | time > 655.715 & time < 656.01;
bad_Ts_dl = bad_Ts_dl | time > 789.72 & time < 790.01;
bad_Ts_dl = bad_Ts_dl | time > 1074.95 & time < 1075.015;
bad_Ts_dl = bad_Ts_dl | time > 1171.76 & time < 1172.01;
bad_Ts_dl = bad_Ts_dl | time > 1288.13 & time < 1289.01;
bad_Ts_dl = bad_Ts_dl | time > 1414.03 & time < 1415.01;
bad_Ts_dl = bad_Ts_dl | time > 1550.88 & time < 1551.01;
bad_Ts_dl = bad_Ts_dl | time > 1580.894 & time < 1580.898;
bad_Ts_dl = bad_Ts_dl | time > 1585.88 & time < 1586.005;
bad_Ts_dl = bad_Ts_dl | time > 1749.9 & time < 1750.005;
bad_Ts_dl = bad_Ts_dl | time > 1838.05 & time < 1839.01;
bad_Ts_dl = bad_Ts_dl | time > 1861.05 & time < 1861.505;
D(84,bad_Ts_dl)=NaN;

%clean GOES
Ts_goes_isnan=isnan(D(210,:))==1 | isempty(D(210,:)) ==1 | isfinite(D(210,:)) == 0;
bad_Ts_goes = Ts_goes_isnan ==1;
bad_Ts_goes = bad_Ts_goes | time > 699.85 & time < 699.9;
bad_Ts_goes = bad_Ts_goes | time > 1082.6 & time < 1083.7;
bad_Ts_goes = bad_Ts_goes | time > 1193.874 & time < 1193.876;
bad_Ts_goes = bad_Ts_goes | time > 1254.68 & time < 1254.7;
bad_Ts_goes = bad_Ts_goes | time > 1294.68 & time < 1294.7;
bad_Ts_goes = bad_Ts_goes | time > 1368.956 & time < 1368.96;
focus = time > 1428.5 & time < 1438.15;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 1441.5 & time < 1700;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 1700 & time < 2150;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 2100 & time < 2220;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 2238 & time < 2238.3;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
D(210,bad_Ts_goes)=NaN;
focus = time > 2290 & time < 2430;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 2515 & time < 2542;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 2581 & time < 2584;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 2620 & time < 2800;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
focus = time > 2564 & time < 2565;
bad_Ts_goes = bad_Ts_goes | focus & D(210,:) == 0;
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
bad = time > 132.67 & time < 133.01;
bad = bad | time > 174.62 & time < 174.845;
bad = bad | time > 220.9 & time < 220.96;
bad = bad | time > 332.7 & time < 333.01;
bad = bad | time > 382.965 & time < 383.01;
bad = bad | time > 554.3 & time < 555.01;
bad = bad | time > 655.72 & time < 655.84;
bad = bad | time > 789.71 & time < 790.01;
bad = bad | time > 1074.95 & time < 1075.01;
bad = bad | time > 1171.76 & time < 1172.01;
bad = bad | time > 1414.035 & time < 1415.01;
bad = bad | time > 1585.885 & time < 1586.01;
bad = bad | time > 1746.7 & time < 1747.005;
bad = bad | time > 1861.05 & time < 1861.51;
bad = bad | time > 1887.95 & time < 1888.01;
focus = time > 2172 & time < 2182;
bad = bad | focus & D(98,:) < 0.02;
D(98,bad) = NaN;

%goes
bad = time > 554.51 & time < 554.97;
bad = bad | time > 1181.075 & time < 1181.09;
bad = bad | time > 1193.8745 & time < 1193.8751;
bad = bad | time > 1254.68 & time < 1254.7;
focus = time > 2170& time < 2180;
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
bad_H2Oc_mat = bad_H2Oc_mat |D(42,:) > 29.8;
bad_H2Oc_mat = bad_H2Oc_mat | time > 1819.665 & time < 1819.67;
focus = time > 2560 & time < 2690;
bad_H2Oc_mat = bad_H2Oc_mat | (focus & D(42,:) > 24);
D(42, bad_H2Oc_mat) = NaN;

%clean data logger
H2Oc_dl_isnan = isnan(D(86,:))==1 | isempty(D(86,:)) ==1 | isfinite(D(86,:)) == 0;
bad_H2Oc_dl = H2Oc_dl_isnan ==1;
bad_H2Oc_dl = bad_H2Oc_dl | D(86,:) < 0;
bad_H2Oc_dl = bad_H2Oc_dl | D(86,:) > 29.8;
bad_H2Oc_dl = bad_H2Oc_dl | time > 801.88 & time < 801.95;
bad_H2Oc_dl = bad_H2Oc_dl | time > 1453 & time < 1462;
bad_H2Oc_dl = bad_H2Oc_dl | time > 1503 & time < 1529;
focus = time > 2560 & time < 2690;
bad_H2Oc_dl = bad_H2Oc_dl | (focus & D(86,:) > 24);
D(86, bad_H2Oc_dl) = NaN;

%convert to mixing ratio
mixing_dl=D(86, :)./(1-(D(86, :)/1000)); %mol h2o/mol dry air
D(86, :)= mixing_dl; %mmol h2o/mol dry air

%clean GOES
%clean GOES
H2Oc_goes_isnan = isnan(D(212,:))==1 | isempty(D(212,:)) ==1 | isfinite(D(212,:)) == 0;
bad_H2Oc_goes = H2Oc_goes_isnan ==1;
bad_H2Oc_goes = bad_H2Oc_goes | D(212,:) < 0;
bad_H2Oc_goes = bad_H2Oc_goes | D(212,:) > 29.8; % The measurement saturated at 30 mmol h2o/ mol
bad_H2Oc_goes = bad_H2Oc_goes | time > 1453 & time < 1462;
bad_H2Oc_goes = bad_H2Oc_goes | time > 1453 & time < 1462;
bad_H2Oc_goes = bad_H2Oc_goes | time > 1284.39 & time < 1284.4;
focus = time > 800 & time < 810;
bad_H2Oc_goes = bad_H2Oc_goes | focus & D(212,:) == 0;
focus = time > 2600 & time < 2700;
bad_H2Oc_goes = bad_H2Oc_goes | (focus & D(212,:) > 24);
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
Tkact_mat   = (D(8,:)+273.15)  ./ (1 + 0.00032 * h2o_mol_frac); %[C]
Tkact_dl   = (D(84,:)+273.15)  ./ (1 + 0.00032 * h2o_mol_frac); %[C]
Tkact_goes   = (D(210,:)+273.15)  ./ (1 + 0.00032 * h2o_mol_frac); %[C]

Tkact_mat   =  Tkact_mat - 273.15; %[C]
Tkact_dl   =  Tkact_dl - 273.15; %[C]
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
%~10% of observations willbe corrected
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
bad = bad | D(21,:) > 1.6; 
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

bad = D(29,:)<-200;
bad = bad | D(29,:)> 550;
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
bad_fastCO2 = D(36,:) == 0; %co2 variability = 0
bad_fastCO2 = bad_fastCO2 | isnan(D(36,:)) == 1; %missing co2 measurements
bad_fastCO2 = bad_fastCO2 | D(35,:) < 340; % remove low concentrations
bad_fastCO2 = bad_fastCO2 | time > 132.6 & time < 133.2;
bad_fastCO2 = bad_fastCO2 | time > 174.682 & time < 174.8;
bad_fastCO2 = bad_fastCO2 | time > 207.8 & time < 208.01;
bad_fastCO2 = bad_fastCO2 | time > 220.92 & time < 220.98;
bad_fastCO2 = bad_fastCO2 | time > 282.975 & time < 283.005;
bad_fastCO2 = bad_fastCO2 | time > 332.8 & time < 333.05;
bad_fastCO2 = bad_fastCO2 | time > 801.92 & time < 801.96;
bad_fastCO2 = bad_fastCO2 | time > 1576.8 & time < 1580.9;
bad_fastCO2 = bad_fastCO2 | time > 1979 & time < 1986.65;
bad_fastCO2 = bad_fastCO2 | time > 1992 & time < 2050;
focus = time > 650 & time < 670;
bad_fastCO2 = bad_fastCO2 | focus & D(36,:) > 6;
D(35,bad_fastCO2)=NaN;

%clean data logger
bad_dlCO2 = D(85,:) < 340; % remove low concentrations
bad_dlCO2 = bad_dlCO2 | time > 132.6 & time < 133.2;
bad_dlCO2 = bad_dlCO2 | time > 174.6 & time < 174.8;
bad_dlCO2 = bad_dlCO2 | time > 207.8 & time < 208.01;
bad_dlCO2 = bad_dlCO2 | time > 220.9 & time < 221;
bad_dlCO2 = bad_dlCO2 | time > 282.88 & time < 282.98;
bad_dlCO2 = bad_dlCO2 | time > 332.7 & time < 333.05;
bad_dlCO2 = bad_dlCO2 | time > 1288 & time < 1289;
bad_dlCO2 = bad_dlCO2 | time > 1414.04 & time < 1415.01;
bad_dlCO2 = bad_dlCO2 | time > 1576.8 & time < 1580.9;
bad_dlCO2 = bad_dlCO2 | time > 1749.93 & time < 1750.01;
bad_dlCO2 = bad_dlCO2 | time > 1979 & time < 1986.65;
bad_dlCO2 = bad_dlCO2 | time > 1992 & time < 2050;
focus = time > 1400 & time < 2700;
bad_dlCO2 = bad_dlCO2 | focus & D(85,:) <385;
D(85, bad_dlCO2) = NaN;

D(85,:)=D(85,:)./(1-(h2o_mol_frac(1,:)/1000)); %co2/mol dry air

%clean GOES
bad_goesCO2 = D(211,:) < 340; % remove low concentrations
bad_goesCO2 = bad_goesCO2 | time > 1576.8 & time < 1580.9;
bad_goesCO2 = bad_goesCO2 | time > 1677.12 & time < 1677.135;
bad_goesCO2 = bad_goesCO2 | time > 1928.1 & time < 1935.98;
bad_goesCO2 = bad_goesCO2 | time > 1966 & time < 2047.98;
focus = time > 1400 & time < 2770;
bad_goesCO2 = bad_goesCO2 | focus & D(211,:) <385;
D(211, bad_goesCO2) = NaN;

D(211,:)=D(211,:)./(1-(h2o_mol_frac(1,:)/1000)); %co2/mol dry air
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
bad_fco2_dl = bad_fco2_dl | time > 659.3 & time < 662; % sensor not working
D(63, bad_fco2_dl) = NaN;


%clean GOES
bad_fco2_goes = bad_goesCO2;
bad_fco2_goes = bad_fco2_goes | isnan(D(186,:))==1; % sonic is not working no uz uz
bad_fco2_goes = bad_fco2_goes | D(186,:)==0;%no wind variability
bad_fco2_goes = bad_fco2_goes | D(201,:)==0; % co2 variability
bad_fco2_goes = bad_fco2_goes | isnan(D(201,:))==1; % sonic is not working no uz uz
bad_fco2_goes = bad_fco2_goes | isnan(D(189,:)) == 1; % no data
bad_fco2_goes = bad_fco2_goes | time > 659.3 & time < 662; % sensor not working
D(189, bad_fco2_goes) = NaN;

%----------------------------------------------------------------------
%Adjust water fluxes using IRGA H2O variabilty compared with HMP variability
%NOTE: IRGA H2O variabilty seems high compared with HMP variability.
%Therefore, we over estimate FH2O.  This will reduce the water flux.

%Make sure the water vapor concentrations are cleaned first.
%----------------------------------------------------------------------
%IRGA calibrations
%The following is the data from the SiteDetails xls spreadsheet.  - Checked on 11/27/2012 awf
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
    outrageous_p = isnan(p(1,1))==1 | isfinite(p(1,1)) == 0 | p(1,1) > 1.3 | p(1,1) < 0.8
    
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
    
    outrageous_p = isnan(p(1,1))==1 | isfinite(p(1,1)) == 0 | p(1,1) > 1.3 | p(1,1) < 0.8
    
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
    outrageous_p = isnan(p(1,1))==1 | isfinite(p(1,1)) == 0 | p(1,1) > 1.3 | p(1,1) < 0.8
    
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
save('..\combined\CombinedInfo\DC_Burn\HMP_h2o_calibration_DC_Burn.mat', 'Pmat', 'Pdl', 'Pgoes');
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
bad = time > 245.84 & time < 254; % 
bad = bad | time > 332 & time < 332.8;
bad = bad | time > 343 & time < 376.7; % bad Rn dome
bad = bad | time > 708.5 & time < 709.5; 
bad = bad | time > 806.5 & time < 807.5; 
%bad = bad | time > 700 & time < 964.92; % Rn may have been bad since ~ day 700
bad = bad | time > 879 & time < 964.92; % bad Rn dome
bad = bad | time > 1436 & time < 1537.86; % bad Rn dome
bad = bad | time > 2582 & time < 2630; % bad Rn dome
D(87,bad) = NaN;

%clean GOES
bad = time > 245.8 & time < 254; % 
bad = bad | time > 343 & time < 376.7; % bad Rn dome
bad = bad | time > 708.5 & time < 709.5; 
bad = bad | time > 806.5 & time < 807.5; 
bad = bad | time > 879 & time < 964.92; % bad Rn dome
bad = bad | time > 1436 & time < 1537.86; % bad Rn dome
bad = bad | time > 2582 & time < 2630; % bad Rn dome
D(213,bad) = NaN;

%data logger net radiation shows some problems
%appears to be filling in ingoing with outgoing for some periods
diff=D(213,:)-D(87, :);
bad = diff > 0.5;
bad = bad |  diff < -0.5;
D(87,bad) = NaN;
%----------------------------------------------------------------------
%Incoming solar radiation - pyranometer = SOLAR_IN
%----------------------------------------------------------------------
%clean data logger
bad = time > 380 & time < 382.6;
bad = bad | time > 501.205 & time < 501.215; %
bad = bad | time > 518.65 & time < 552.035; %
bad = bad | time > 719.06 & time < 734.94; %
bad = bad | time > 989.38 & time < 1093.07;
bad = bad | time > 1153.88 & time < 1161.58;
bad = bad | time > 1181.36 & time < 1161.58;
bad = bad | time > 1203.79 & time < 1265;
D(88,bad) = NaN;

%clean goes
bad = time > 380 & time < 383.7;
bad = bad | time > 500.88 & time < 501.215; % 
bad = bad | time > 518.65 & time < 552.035; % 
bad = bad | time > 719.06 & time < 734.94; %
bad = bad | time > 989.38 & time < 1093.07;
bad = bad | time > 1153.88 & time < 1161.58;
bad = bad | time > 1181 & time < 1265.2;
D(214,bad) = NaN;

%data logger radiation shows some problems
diff=D(214,:)-D(88, :);
bad = diff > 1;
bad = bad |  diff < -1;
D(88,bad) = NaN;
%----------------------------------------------------------------------
%Outgoing solar radiation - pyranometer = SOLAR_OUT
%----------------------------------------------------------------------
%clean data logger
bad = D(89,:) < -3;
bad = bad | time > 149.8 & time < 150.84;
bad = bad | time > 198.99 & time < 207.6;
bad = bad | time > 222.005 & time < 242.05;
bad = bad | time > 661.95 & time < 719.05;
bad = bad | time > 809.98 & time < 966.64;
bad = bad | time > 1074.94 & time < 1074.99;
bad = bad | time > 1153.88 & time < 1264.95;
bad = bad | time > 1539.9 & time < 1560.88;
bad = bad | time > 174.7 & time < 174.85;
D(89,bad) = NaN;

%clean goes
bad = D(215,:) < -3;
bad = bad | time > 149.8 & time < 150.84;
bad = bad | time > 174.854 & time < 174.8545;
bad = bad | time > 198.99 & time < 207.6;
bad = bad | time > 220.97 & time < 221.01;
bad = bad | time > 222.005 & time < 227.98;
bad = bad | time > 228.4 & time < 242.05;
bad = bad | time > 661.95 & time < 719.05;
bad = bad | time > 809.98 & time < 966.64;
bad = bad | time > 1153.88 & time < 1264.95;
bad = bad | time > 1539.9 & time < 1560.88;
focus = time > 1900 & time < 1910;
bad = bad | focus & D(215,:) == 0;
D(215,bad) = NaN;

%data logger radiation shows some problems
diff=D(215,:)-D(89,:);
bad = diff > 0.6;
bad = bad |  diff < -0.6;
D(89,bad) = NaN;
%----------------------------------------------------------------------
%PAR_In
%----------------------------------------------------------------------
%clean data logger
bad = D(90,:) < -20; 
bad = bad | time > 131.73 & time < 180;
bad = bad | time > 194.9 & time < 207.56;
bad = bad | time > 208.15 & time < 221.2;
bad = bad | time > 208.15 & time < 221.2;
bad = bad | time > 789.2 & time < 789.72;
bad = bad | time > 790.005 & time < 792.86;
bad = bad | time > 793.86 & time < 808.84;
bad = bad | time > 950.5 & time < 964.5;
bad = bad | time > 1980 & time < 2047.96;
bad = bad | time > 2212.5 & time < 2225.3;
bad = bad | time > 2357 & time < 2571.5;
bad = bad | time > 2584 & time < 2767;
D(90,bad) = NaN;

%clean GOES
bad = D(216,:) < -20;
bad = bad | time > 131.73 & time < 180;
bad = bad | time > 194.9 & time < 207.56;
bad = bad | time > 208.15 & time < 221.2;
bad = bad | time > 789.2 & time < 792.86;
bad = bad | time > 793.86 & time < 808.84;
bad = bad | time > 950.5 & time < 964.5;
bad = bad | time > 1980 & time < 2047.96;
bad = bad | time > 2212.19 & time < 2225.5;
bad = bad | time > 2356.06 & time < 2571.5;
bad = bad | time > 2584 & time < 2767;
D(216,bad) = NaN;

%----------------------------------------------------------------------
%PAR_Out
%----------------------------------------------------------------------
%use this for cleaning night time 
hour = time-floor(time);

%clean data logger
bad = D(91,:) <-8; 
bad = bad | D(91,:) > 1000;
bad = bad | time > 174.7 & time < 174.85;
bad = bad | time > 207 & time < 221.03;
bad = bad | time > 317.6 & time < 324.5;
bad = bad | time > 325.87 & time < 325.88;
bad = bad | time > 377.6 & time < 378;
bad = bad | time > 452.7 & time < 452.9;
bad = bad | time > 467.3 & time < 470.45;
bad = bad | time > 719.5 & time < 722.8;
bad = bad | time > 776.6 & time < 776.8;
bad = bad | time > 870.74 & time < 941.01;
bad = bad | time > 998.91 & time < 1153.94;
bad = bad | time > 1155.78 & time < 1165.7;
bad = bad | time > 1167.75 & time < 1167.8;
bad = bad | time > 1169.3 & time < 1169.9;
bad = bad | time > 1172.815 & time < 1173.28;
bad = bad | time > 1173.82 & time < 1182.4;
bad = bad | time > 1183.93 & time < 1184.1;
bad = bad | time > 1190 & time < 1265;
focus = time > 1307.4 & time < 1533;
bad = bad | focus & D(91,:) > 320;
bad = bad | time > 1309.1 & time < 1311.24;
bad = bad | time > 1335 & time < 1335.4;
bad = bad | time > 1335.51 & time < 1335.53;
bad = bad | time > 1395.4 & time < 1395.6;
bad = bad | time > 1408.1 & time < 1410;
bad = bad | time > 1420.35 & time < 1420.36;
bad = bad | time > 1420.88 & time < 1421.96;
bad = bad | time > 1433 & time < 1440;
bad = bad | time > 1529 & time < 1531;
focus = time > 1538.2 & time < 1538.65;
bad = bad | focus & D(91,:) > 150;
focus = time > 1538.2 & time < 1551;
focus = focus & hour > 0.12 & hour < 0.55;
bad = bad | focus & D(91,:) > 2;
bad = bad | time > 1553.71 & time < 1553.84;
bad = bad | time > 1559.5 & time < 1560;
D(91,bad) = NaN;

%clean GOES
bad = D(217,:) <-8; 
bad = bad | D(217,:) > 1000;
bad = bad | time > 147.93 & time < 148.06;
bad = bad | time > 174.85 & time < 174.86;
bad = bad | time > 207 & time < 221.03;
bad = bad | time > 227.96 & time < 228.14;
bad = bad | time > 317.6 & time < 324.5;
focus = time > 325.86 & time < 330;
bad = bad | focus & D(217,:) == 0;
bad = bad | time > 325.87 & time < 325.88;
bad = bad | time > 374.85 & time < 375.05;
bad = bad | time > 375.76 & time < 376;
bad = bad | time > 377.6 & time < 378;
bad = bad | time > 452.7 & time < 452.9;
bad = bad | time > 467.3 & time < 470.45;
bad = bad | time > 492.9 & time < 492.96;
bad = bad | time > 719.5 & time < 722.8;
bad = bad | time > 776.6 & time < 776.8;
bad = bad | time > 870.74 & time < 941.01;
bad = bad | time > 998.91 & time < 1153.94;
bad = bad | time > 1155.78 & time < 1165.7;
bad = bad | time > 1167.75 & time < 1167.8;
bad = bad | time > 1169.3 & time < 1169.9;
bad = bad | time > 1172.815 & time < 1173.28;
bad = bad | time > 1173.82 & time < 1182.4;
bad = bad | time > 1183.93 & time < 1184.1;
focus = time > 1186 & time < 1190;
bad = bad | focus & D(217,:) == 0;
focus = time > 1303 & time < 1305;
bad = bad | focus & D(217,:) == 0;
bad = bad | time > 1190 & time < 1265;
bad = bad | time > 1308.9 & time < 1310.46;
bad = bad | time > 1310.72 & time < 1310.74;
bad = bad | time > 1311.15 & time < 1311.25;
bad = bad | time > 1318.6 & time < 1319;
bad = bad | time > 1335 & time < 1335.4;
bad = bad | time > 1337.5 & time < 1337.7;
bad = bad | time > 1395.4 & time < 1395.6;
bad = bad | time > 1408.3 & time < 1409.7;
bad = bad | time > 1420.35 & time < 1420.36;
bad = bad | time > 1420.88 & time < 1421.96;
bad = bad | time > 1433 & time < 1440;
bad = bad | time > 1529 & time < 1531;
bad = bad | time > 1553.71 & time < 1553.84;
focus = time > 1307.4 & time < 1533;
bad = bad | focus & D(217,:) > 320;
bad = bad | focus & D(217,:) == 0;
focus = time > 1534.5 & time < 1561;
bad = bad | focus & D(217,:) == 0;
bad = bad | time > 1538.3 & time < 1538.635;
bad = bad | time > 1559.2 & time < 1560.4;
focus = time > 1538.2 & time < 1551;
focus = focus & hour > 0.12 & hour < 0.55;
bad = bad | focus & D(91,:) > 2;
bad = bad | time > 1877.6 & time < 1878;
D(217,bad) = NaN;

%data logger radiation shows some problems
diff=D(217,:)-D(91,:);
bad = diff > 5;
bad = bad |  diff < -5;
D(91,bad) = NaN;
%----------------------------------------------------------------------
%T_107
%----------------------------------------------------------------------
%dl
D(102,:) = NaN;

%goes
D(228,:) = NaN;

%----------------------------------------------------------------------
%Rain
%----------------------------------------------------------------------
%data logger
bad = D(103,:)>40;
D(103,bad) = NaN;

%goes
bad = D(229,:)>40;
D(229,bad) = NaN;

%----------------------------------------------------------------------
%NDVI
%----------------------------------------------------------------------

%----------------------------------------------------------------------
%Soil T
%----------------------------------------------------------------------

%T1
bad = D(288,:) < -10;
D(288,bad)=NaN;

%T2
bad = D(289,:) < -10;
D(289,bad)=NaN;

%T3
bad = D(290,:) < -10;
bad = bad | time > 1804.7 & time < 1893.76;
D(290,bad)=NaN;

%T4
bad = D(291,:) < -10;
D(291,bad)=NaN;

%----------------------------------------------------------------------
%Clean Soil Moisture
%----------------------------------------------------------------------
%looks good - awf

%----------------------------------------------------------------------
%Fuel_M 
%----------------------------------------------------------------------
bad = D(286,:) < 4;
D(286,bad)=NaN;
%----------------------------------------------------------------------
%LWS
%----------------------------------------------------------------------
%don't really know what this data should look like
bad = time > 1685.17 & time < 2695;
D(295,bad)=NaN;

bad = time > 2254.215 & time < 2695;
D(294,bad)=NaN;

bad = D(293,:) < 220;
D(293,bad)=NaN;
%----------------------------------------------------------------------
%Clean Matric potential sensors 
%----------------------------------------------------------------------
bad = time > 1331.06 & time < 1331.1;
bad = bad | D(314,:) > 3;
bad = bad | D(314,:) < 0;
D(314,bad)=NaN;

bad = time > 1331.06 & time < 1331.1;
bad = bad | D(315,:) > 3;
bad = bad | D(315,:) < 0;
D(315,bad)=NaN;

bad = time > 1331.06 & time < 1331.1;
bad = bad | D(316,:) > 3;
bad = bad | D(316,:) < 0;
D(316,bad)=NaN;

bad = time > 1331.06 & time < 1331.1;
bad = bad | D(317,:) > 3;
bad = bad | D(317,:) < 0;
D(317,bad)=NaN;

bad = time > 1331.06 & time < 1331.1;
bad = bad | D(318,:) > 3;
bad = bad | D(318,:) < 0;
D(318,bad)=NaN;

bad = time > 1331.06 & time < 1331.1;
bad = bad | D(319,:) > 3;
bad = bad | D(319,:) < 0;
D(319,bad)=NaN;

elseif a==2
  
%Mop-up

TEMP =[];
H2O =[];
rho_fix =[];
v_wind = [];    
  
%Sensible Heat
%dl
bad = D(339,:)<-200;
bad = bad | D(339,:)> 550;
D(339,bad)=NaN;

%goes
bad = D(340,:)<-200;
bad = bad | D(340,:)> 550;
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
bad_le_dl = bad_dl_et;
D(343,bad_le_dl)=NaN;

%clean GOES 
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