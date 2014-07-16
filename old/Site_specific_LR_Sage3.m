function [D] = Site_specific_LR_Sage(HEADER, D)

%This is a site specific file for Loma Ridge Grassland

%use time to remove bad data
time = D(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%By measurement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------
%Tsonic
%----------------------------------------------------------------------
%clean mat
bad_Ts_mat=D(8,:)<1;
D(8,bad_Ts_mat)=NaN;

%clean dl
bad_Ts_dl=D(84,:)<1;
D(84,bad_Ts_dl)=NaN;

%clean GOES
bad_Ts_goes=D(210,:)<1;
D(210,bad_Ts_goes)=NaN;

%----------------------------------------------------------------------
%Wind Speed
%----------------------------------------------------------------------
%If wind vectors are wrong, then all wind and momentum fluxes are wrong
%clean mat
%component wind problems
u_mat=D(5,:);
v_mat=D(6,:);
w_mat=D(7,:);
bad_wind_mat=u_mat==0;
bad_wind_mat = bad_wind_mat | v_mat==0;
bad_wind_mat = bad_wind_mat | w_mat==0;

%clean dl
%DATA LOGGER
u_dl=D(81,:);%'Ux_1_Avg(1)');
v_dl=D(82,:);%'Uy_1_Avg(1)');
w_dl=D(83,:);
bad_wind_dl=u_dl==0;
bad_wind_dl = bad_wind_dl | v_dl==0;
bad_wind_dl = bad_wind_dl | w_dl==0;

%clean goes
u_goes = D(207,:);%'GOES_U_avg'); 
v_goes = D(208,:);%'GOES_Uy_avg'); 
w_goes = D(209,:);
bad_wind_goes = u_goes == 0;
bad_wind_goes = bad_wind_goes | v_goes==0;
bad_wind_goes = bad_wind_goes | w_goes==0;


%clean fast
D(327,bad_wind_mat)=NaN;

%clean data logger
D(328,bad_wind_dl)=NaN;

%clean GOES
D(329,bad_wind_goes) = NaN; 


%----------------------------------------------------------------------
%Ustar
%----------------------------------------------------------------------
%clean FAST
bad = bad_wind_mat;
bad = bad | D(21,:) > 1.5;
D(21,bad)=NaN;

%clean data logger
bad = bad_wind_dl;
bad = bad | D(61,:) == 0; %covariance Uz_Ux
bad = bad | D(62,:) == 0; %covariance Uz_Ux
bad = bad | D(322,:) > 1.5;
D(322,bad)=NaN;

%clean GOES
bad = bad_wind_goes;
bad = bad | D(187,:) == 0; %covariance Uz_Ux
bad = bad | D(188,:) == 0; %covariance Uz_Ux
bad = bad | D(323,:) > 1.5;
D(323,bad)=NaN;

%----------------------------------------------------------------------
%Wind Direction
%----------------------------------------------------------------------
D(324,bad_wind_mat) = NaN; %fast
D(325,bad_wind_dl) = NaN;%dl
D(326,bad_wind_goes) = NaN;%goes

%----------------------------------------------------------------------
%bad sonic
%----------------------------------------------------------------------
bad_uvwt_mat=bad_wind_mat|bad_Ts_mat;
bad_uvwt_dl=bad_wind_dl|bad_Ts_dl;
bad_uvwt_goes=bad_wind_goes|bad_Ts_goes;


%----------------------------------------------------------------------
%Tactual
%----------------------------------------------------------------------
%clean mat
D(330,bad_Ts_mat)=NaN;

%clean dl
D(331,bad_Ts_dl)=NaN;

%clean GOES
D(332,bad_Ts_goes)=NaN;

%----------------------------------------------------------------------
%Dry Density
%----------------------------------------------------------------------
%If dry density is wrong then the fluxes are wrong too
%clean mat
D(336,bad_Ts_mat)=NaN;

%clean dl
D(337,bad_Ts_dl)=NaN;

%clean GOES
D(338,bad_Ts_goes)=NaN;

%----------------------------------------------------------------------
%Moist Density
%----------------------------------------------------------------------
%clean mat
D(333,bad_Ts_mat)=NaN;

%clean dl
D(334,bad_Ts_dl)=NaN;

%clean GOES
D(335,bad_Ts_goes)=NaN;

%----------------------------------------------------------------------
%Sensible heat
%----------------------------------------------------------------------
%clean fast
bad = bad_uvwt_mat;
bad = bad | D(29,:)<-300;
D(29,bad)=NaN;

%clean data logger
bad = bad_uvwt_dl;
bad = bad | D(339,:)<-300;
bad = bad | D(339,:)> 1100;
bad = bad | time > 1812.9374 & time < 1812.9376;
D(339,bad)=NaN;

%clean GOES
bad = bad_uvwt_goes;
bad = bad | D(340,:)<-300;
bad = bad | D(340,:)> 1100;
D(340, bad) = NaN;

%----------------------------------------------------------------------
%CO2 concentration
%----------------------------------------------------------------------
%clean fast
bad_fastCO2 = D(35,:) < 360;
D(35,bad_fastCO2)=NaN;

%clean data logger
bad_dlCO2 = D(85,:) < 360;
D(85, bad_dlCO2) = NaN;

%clean GOES
bad_goesCO2 = D(211,:) < 360;
D(211, bad_goesCO2) = NaN;

%----------------------------------------------------------------------
%FCO2
%----------------------------------------------------------------------
%clean fast
D(46, bad_fastCO2) = NaN;
D(46, bad_uvwt_mat) = NaN;

%out of reasonable range range
bad = D(46,:) < -50;
bad = bad | D(46,:) > 30;
bad = bad | D(46,:) == 0;
D(46, bad) = NaN;

%clean data logger
D(345, bad_dlCO2) = NaN;
D(345, bad_uvwt_dl) = NaN;

bad = D(345,:) < -50;
bad = bad | D(345,:) > 30;
bad = bad | D(345,:) == 0;
D(345, bad) = NaN;

%clean GOES
D(346, bad_goesCO2) = NaN;
D(346, bad_uvwt_goes) = NaN;

bad = D(346,:) < -50;
bad = bad | D(346,:) > 30;
bad = bad | D(346,:) == 0;
D(346, bad) = NaN;

%----------------------------------------------------------------------
%IRGA H2O concentrations 
%carry bad IRGA H2O concentrations  through the water fluxes
%----------------------------------------------------------------------
%clean fast
bad_H2Oc_mat = D(42,:) <= 0;
bad_H2Oc_mat = bad_H2Oc_mat | time > 1080.583 & time < 1080.7917;
D(42, bad_H2Oc_mat) = NaN;

%clean data logger
bad_H2Oc_dl = D(86,:) <= 0;
D(86, bad_H2Oc_dl) = NaN;

%clean GOES
bad_H2Oc_goes = D(212,:) <= 0;
D(212, bad_H2Oc_goes) = NaN;


%----------------------------------------------------------------------
%FH2O
%----------------------------------------------------------------------
%clean fast
D(47,bad_H2Oc_mat)=NaN;
D(47, bad_uvwt_mat) = NaN;

bad = bad | D(47,:) < -10;
bad = bad | D(47,:) > 15;
bad = bad | D(47,:) == 0;
D(47, bad) = NaN;

%clean data logger
%clean data logger
D(347, bad_H2Oc_dl) = NaN;
D(347, bad_uvwt_dl) = NaN;

bad = bad | D(347,:) < -10;
bad = bad | D(347,:) > 15;
bad = bad | D(347,:) == 0;
D(347, bad) = NaN;

%clean GOES
D(348, bad_H2Oc_goes) = NaN;
D(348, bad_uvwt_goes) = NaN;

bad = bad | D(348,:) < -10;
bad = bad | D(348,:) > 15;
bad = bad | D(348,:) == 0;
D(348, bad) = NaN;

%----------------------------------------------------------------------
%Latent heat
%----------------------------------------------------------------------
%fast
D(30,bad_H2Oc_mat)=NaN;
D(30, bad_uvwt_mat) = NaN;

bad = D(30,:) <-300;
bad= bad | D(30,:) >500;
bad= bad | D(30,:) == 0;
D(30,bad)=NaN;

%clean data logger
D(343, bad_H2Oc_dl) = NaN;
D(343, bad_uvwt_dl) = NaN;

bad = D(343,:) <-300;
bad= bad | D(343,:) >500;
bad= bad | D(343,:) == 0;
D(343,bad)=NaN;


%clean GOES
D(344, bad_H2Oc_goes) = NaN;
D(344, bad_uvwt_goes) = NaN;

bad = D(344,:) <-300;
bad= bad | D(344,:) >500;
bad= bad | D(344,:) == 0;
D(344,bad)=NaN;

%----------------------------------------------------------------------
%Rn
%----------------------------------------------------------------------
%clean data logger
bad = D(87,:) == 0; 
bad = bad | time > 545.8543 & time < 545.938; 
D(87,bad) = NaN;

%clean GOES
bad = D(213,:) == 0; 
bad = bad | D(213,:) > 2000; 
bad = bad | time < 88.6; 
D(213,bad) = NaN;

%----------------------------------------------------------------------
%Incoming solar radiation - pyranometer = SOLAR_IN
%----------------------------------------------------------------------
%clean data logger
bad = D(88,:) == 0; 
bad = bad | time > 545.8543 & time < 545.938; 
D(88,bad) = NaN;

%clean goes
bad = D(214,:) == 0; 
bad = bad | time > 1178.0416 & time < 1178.0418; 
D(214,bad) = NaN;

%----------------------------------------------------------------------
%Outgoing solar radiation - pyranometer = SOLAR_OUT
%----------------------------------------------------------------------
%clean data logger
bad = D(89,:) == 0; 
bad = bad | D(89,:) < -20; 
bad = bad | D(89,:) > 300; 
bad = bad | time > 980 & time < 1228; 
D(89,bad) = NaN;

%clean goes
bad = D(215,:) == 0; 
bad = bad | D(215,:) < -20; 
bad = bad | D(215,:) > 300; 
bad = bad | time > 980 & time < 1228; 
D(215,bad) = NaN;

%----------------------------------------------------------------------
%PAR_In
%----------------------------------------------------------------------
%clean data logger
bad = D(90,:) == 0; 
bad = bad | time > 545.8543 & time < 545.938; 
bad = bad | time > 1012 & time < 1228; 
D(90,bad) = NaN;

%clean GOES
bad = D(216,:) == 0; 
bad = bad | time > 1012 & time < 1228; 
D(216,bad) = NaN;

%----------------------------------------------------------------------
%PAR_Out
%----------------------------------------------------------------------
%clean data logger
bad = D(91,:) == 0; 
bad = bad | D(91,:) < -20;
bad = bad | time(1,:) > 104 & time(1,:) < 119.9;
bad = bad | time(1,:) > 1260 & time(1,:) < 1355.9995;
D(91,bad) = NaN;

%clean GOES
bad = D(217,:) == 0; 
bad = bad | D(217,:) < -20;
bad = bad | D(217,:) > 1000;
bad = bad | time(1,:) > 104 & time(1,:) < 119.9;
bad = bad | time(1,:) > 1260 & time(1,:) < 1355.9995;
D(217,bad) = NaN;

%----------------------------------------------------------------------
%HMP_Temp 
%----------------------------------------------------------------------
%data logger
bad=D(97,:)<1;
bad= bad | time > 458 & time < 460;
D(97,bad) = NaN;

%goes
bad=D(223,:)<1;
bad= bad | time > 1599.7 & time < 1601.1;
D(223,bad) = NaN;

%----------------------------------------------------------------------
%Water vapor concentration
%----------------------------------------------------------------------
%data logger
bad = time(1,:) > 458.55 & time(1,:) < 518;
D(100,bad) = NaN;

%goes
bad = time(1,:) > 458.55 & time(1,:) < 518;
D(226,bad) = NaN;
%----------------------------------------------------------------------
%HMP_RH
%----------------------------------------------------------------------
%data logger
bad = time(1,:) > 458.55 & time(1,:) < 518;
bad= bad | D(98,:) <= 0;
D(98,bad) = NaN;

%goes
bad = time(1,:) > 458.55 & time(1,:) < 518;
bad= bad | time > 1599.7 & time < 1601.1;
bad= bad | D(224,:) <= 0;
D(224,bad) = NaN;

%----------------------------------------------------------------------
%T_107
%----------------------------------------------------------------------
%dl
bad = D(102,:) < 1;
D(102,bad) = NaN;

%goes
bad = D(228,:) < 1;
bad= bad | time > 1599.7 & time < 1601.1;
D(228,bad) = NaN;


%----------------------------------------------------------------------
%Rain
%----------------------------------------------------------------------
%data logger
bad=D(103,:)>5000;
D(103,bad) = NaN;


D=D;


