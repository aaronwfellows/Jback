function [D] = Site_specific_LR_Grass(HEADER, D)

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
%Tactual
%----------------------------------------------------------------------
%clean GOES
bad_Tactual=D(332,:)<274;
D(332,bad_Tactual)=NaN;

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
D(21,bad)=NaN;

%clean data logger
bad = bad_wind_dl;
bad = bad | D(61,:) == 0; %covariance Uz_Ux
bad = bad | D(62,:) == 0; %covariance Uz_Ux
bad = bad | time > 858.2918 & time < 858.66;
D(322,bad)=NaN;

%clean GOES
bad = bad_wind_goes;
bad = bad | D(187,:) == 0; %covariance Uz_Ux
bad = bad | D(188,:) == 0; %covariance Uz_Ux
bad = bad | time > 858.2918 & time < 858.66;
D(323,bad)=NaN;

%----------------------------------------------------------------------
%Wind Direction
%----------------------------------------------------------------------
D(324,bad_wind_mat) = NaN; %fast
D(325,bad_wind_dl) = NaN;%dl
D(326,bad_wind_goes) = NaN;%goes

%----------------------------------------------------------------------
%Wind Speed
%----------------------------------------------------------------------
%clean fast
D(327,bad_wind_mat)=NaN;

%clean data logger
D(328,bad_wind_dl)=NaN;

%clean GOES
D(329,bad_wind_goes) = NaN; 


%----------------------------------------------------------------------
%bad sonic
%----------------------------------------------------------------------
bad_uvwt_mat=bad_wind_mat|bad_Ts_mat;
bad_uvwt_dl=bad_wind_dl|bad_Ts_dl;
bad_uvwt_goes=bad_wind_goes|bad_Ts_goes;

%----------------------------------------------------------------------
%Dry Density
%----------------------------------------------------------------------
%If dry density is wrong then the fluxes are wrong too
   
%clean GOES
%GOES has problems due to spiking Tactual used to calculate the density
D(338,bad_Tactual)=NaN;

%----------------------------------------------------------------------
%Moist Density
%----------------------------------------------------------------------
%clean GOES
%GOES has problems due to spiking Tactual used to calculate the density
D(335,bad_Tactual)=NaN;

%----------------------------------------------------------------------
%Sensible heat
%----------------------------------------------------------------------
%clean fast
D(29,bad_uvwt_mat)=NaN;
bad=D(29,:)<300;
D(29,bad)=NaN;

%clean data logger
D(339, bad_uvwt_dl) = NaN;

%clean GOES
D(340, bad_uvwt_goes) = NaN;

%----------------------------------------------------------------------
%CO2 concentration
%----------------------------------------------------------------------
%clean fast
bad_fastCO2 = D(35,:) < 344;
bad_fastCO2 = bad_fastCO2 | time > 1537.9165 & time < 1537.917;
D(35,bad_fastCO2)=NaN;

%clean data logger
bad_dlCO2 = D(85,:) < 344;
D(85, bad_dlCO2) = NaN;

%clean GOES
bad_goesCO2 = D(211,:) < 344;
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
bad_H2Oc_mat = D(42,:) < 0.01;
D(42, bad_H2Oc_mat) = NaN;

%clean data logger
bad_H2Oc_dl = D(86,:) < 0.01;
bad_H2Oc_dl = bad_H2Oc_dl | time(1,:) > 510 & time(1,:) < 540 & D(86,:) < 3;
D(86, bad_H2Oc_dl) = NaN;

%clean GOES
bad_H2Oc_goes = D(212,:) < 0.01;
bad_H2Oc_goes = bad_H2Oc_goes | time(1,:) > 510 & time(1,:) < 540 & D(212,:) < 3;
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
D(30,bad)=NaN;
bad= D(30,:) >500;
D(30,bad)=NaN;

%clean data logger
D(343, bad_H2Oc_dl) = NaN;
D(343, bad_uvwt_dl) = NaN;

bad = D(343,:) <-300;
D(343,bad)=NaN;
bad= D(343,:) >500;
D(343,bad)=NaN;

%clean GOES
D(344, bad_H2Oc_goes) = NaN;
D(344, bad_uvwt_goes) = NaN;

bad = D(344,:) <-300;
D(344,bad)=NaN;
bad= D(344,:) >500;
D(344,bad)=NaN;
%----------------------------------------------------------------------
%Rn
%----------------------------------------------------------------------
%clean data logger
bad = D(87,:) == 0; 
bad = bad |time(1,:) > 1706.8 & time(1,:) < 1887.5;
D(87,bad) = NaN;

%clean GOES
bad = D(213,:) == 0; 
bad = bad | time(1,:) > 1706.8 & time(1,:) < 1887.5;
D(213,bad) = NaN;

%----------------------------------------------------------------------
%PAR_In
%----------------------------------------------------------------------
%clean data logger
bad = D(90,:) == 0; 
bad = bad | time(1,:) > 660 & time(1,:) < 892.2;
D(90,bad) = NaN;

%clean GOES
bad = D(216,:) == 0; 
bad = bad | time(1,:) > 660 & time(1,:) < 892.2;
D(216,bad) = NaN;

%----------------------------------------------------------------------
%PAR_Out
%----------------------------------------------------------------------
%clean data logger
bad = D(91,:) == 0; 
bad = bad | time(1,:) > 611.8 & time(1,:) < 706.5;
bad = bad | time(1,:) > 848.97 & time(1,:) < 892.2;
bad = bad | time(1,:) > 1782.6456 & time(1,:) < 1782.646;
D(91,bad) = NaN;

%clean GOES
bad = D(217,:) == 0; 
bad = bad | time(1,:) > 611.8 & time(1,:) < 706.5;
bad = bad | time(1,:) > 848.97 & time(1,:) < 892.2;
D(217,bad) = NaN;

%----------------------------------------------------------------------
%Incoming solar radiation - pyranometer = SOLAR_IN
%----------------------------------------------------------------------
%clean data logger
bad = D(88,:) == 0; 
bad = bad | time(1,:) > 378.74 & time(1,:) < 511;
bad = bad | (time(1,:) > 1358.6 & time(1,:) < 1530.9);
D(88,bad) = NaN;

%clean goes
bad = D(214,:) == 0; 
bad = bad | time(1,:) > 378.74 & time(1,:) < 511;
bad = bad | (time(1,:) > 1358.6 & time(1,:) < 1530.9);
D(214,bad) = NaN;

%----------------------------------------------------------------------
%Outgoing solar radiation - pyranometer = SOLAR_OUT
%----------------------------------------------------------------------

%clean data logger


%clean goes



%----------------------------------------------------------------------
%HMP_Temp 
%----------------------------------------------------------------------
%data logger
bad=D(97,:)<1;
bad = bad | time(1,:) > 587.85 & time(1,:) < 675;
bad = bad | time(1,:) > 1451.8 & time(1,:) < 1531.2;
bad = bad | time(1,:) > 1672.95 & time(1,:) < 1701.5;
D(97,bad) = NaN;

%goes
D(223,bad) = NaN;

bad=D(223,:)<1;
D(223,bad) = NaN;
%----------------------------------------------------------------------
%HMP Water vapor concentration
%----------------------------------------------------------------------
%data logger
bad = time(1,:) > 587.85 & time(1,:) < 675;
bad = bad | time(1,:) > 1451.8 & time(1,:) < 1531.2;
bad = bad | time(1,:) > 1672.95 & time(1,:) < 1701.5;
D(100,bad) = NaN;

%goes
D(226,bad) = NaN;
%----------------------------------------------------------------------
%HMP_RH
%----------------------------------------------------------------------
%data logger
bad = time(1,:) > 587.85 & time(1,:) < 675;
bad = bad | time(1,:) > 1451.8 & time(1,:) < 1531.2;
bad = bad | time(1,:) > 1672.95 & time(1,:) < 1701.5;
D(98,bad) = NaN;

%goes
D(224,bad) = NaN;

%----------------------------------------------------------------------
%T_107
%----------------------------------------------------------------------
%dl
bad = D(102,:) < 1;
bad = bad | time(1,:) < 104.2;
bad = bad | time(1,:) > 1451.8 & time(1,:) < 1531.2;
bad = bad | time(1,:) > 1672.95 & time(1,:) < 1701.5;
D(102,bad) = NaN;

%goes
bad = D(228,:) < 1;
bad = bad | time(1,:) < 104.2;
bad = bad | time(1,:) > 1451.8 & time(1,:) < 1531.2;
bad = bad | time(1,:) > 1672.95 & time(1,:) < 1701.5;
D(228,bad) = NaN;

%----------------------------------------------------------------------
%Rain
%----------------------------------------------------------------------
%data logger
bad=D(103,:)>5000;
D(103,bad) = NaN;

%----------------------------------------------------------------------
%NDVI
%----------------------------------------------------------------------
bad = D(96,:) > 1;
D(96,bad) = NaN;

bad = D(222,:) > 1;
D(222,bad) = NaN;

%----------------------------------------------------------------------
%Soil Temp
%----------------------------------------------------------------------
bad = time > 850.935 & time < 850.96;
D(291,bad) = NaN;
D(292,bad) = NaN;
bad = time > 850.93 & time < 850.95;
D(293,bad) = NaN;
%----------------------------------------------------------------------
%Clean Soil Moisture
%----------------------------------------------------------------------
bad = D(285,:) < 0;
D(285,bad) = NaN;


%----------------------------------------------------------------------
%Clean Matric potential sensors 
%----------------------------------------------------------------------
bad = time(1,:) > 1950 & time(1,:) < 2200;
D(299,bad) = NaN;%T
D(305,bad) = NaN;%T1
D(311,bad) = NaN;%T30
D(317,bad) = NaN;%delT

bad = time(1,:) > 1410 & time(1,:) < 1535.92;
bad = bad | time(1,:) > 2119.41 & time(1,:) < 2200;
D(301,bad) = NaN;
D(307,bad) = NaN;
D(313,bad) = NaN;
D(319,bad) = NaN;

bad = time(1,:) > 1408 & time(1,:) < 1535.92;
D(302,bad) = NaN;
D(308,bad) = NaN;
D(314,bad) = NaN;
D(320,bad) = NaN;

bad = time(1,:) > 1497 & time(1,:) < 2200;
D(309,bad) = NaN;
D(315,bad) = NaN;
D(321,bad) = NaN;

%T1s is bad
bad=D(316,:)<0;
D(316,bad)=NaN;

bad=D(317,:)<0;
D(317,bad)=NaN;

bad=D(318,:)<0;
D(318,bad)=NaN;

%out
D=D;


