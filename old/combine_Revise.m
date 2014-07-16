%function combine_newV3inProgress(sites2Proc)
%What sites do you want to process? see var_defs sites table to relate
%number to site
sites2Proc=4;

%--------------------------------------------------------------------------
%origin:
%Nearly all of this code comes from the Anne Kelly's combined_new program circa
%4/3/12

%Reorginization of code by awf %4/11/12
%--------------------------------------------------------------------------
%The purpose is to combined together 3 eddy covariance data streams: 
%1) fast = processed on computer from 4 Hz data, 2)datalogger = processed on the data logger, 
% and 3)goes = processed on the data logger and transmitted by goes.

%Assumes you processed data using fastflux and put all data
%streams into a big table using the merge program.

%This program:
    % General filtering
    % Sensor calibrations 
    % Each dataset will be filtered individually
    % We then combined these different data streams in the following order: 
    % 1) fast, 2)datalogger, 3)goes to produce 1 data table (Cleaned_D)
    
    %What is the variable? 
    %Cleaned_Header=[];
    
    %Where is the data stored? The data is orginized by column and stored in 
    %Cleaned_D=[];
    
    %What data stream did the datum come from?
    %DataSource=[]; %1 = processed (aka: mat) ; 2=datalogger; 3=GOES
    %DataSource uses the Cleaned_Header

%--------------------------------------------------------------------------
%There are 6 sections
%(1)non-site specific filtering
%(2)Sensor Calibration --> turn to physical units
%(3)Calculations
%(4)Site specific filtering - delves into additional site specific programs
%(5)combined data streams
%(6)save

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Initialilizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Site-specific info
global sites sitePress towerYearStart iSite SonicOrientation VIsites minFlow
% Paths etc.
path(path, 'C:\towerData\processingScripts\subroutines');
% flux processing parameters
global oneAtmc R_mol Tc CpAirc CpH2Oc Mw_da Mw_water
%
    
var_defs(); %this has important constants and tables ***make sure this is up to date***
Day = date;
%%
close all
tic

%diary_filename = [combineRootDir, 'combine_log_', Day, char(sites(sites2Proc(1)))];
%diary(diary_filename);
format long g

%directory to grab merged files and dump combined files
mergedRootDir='C:\towerData\merged\';
combineRootDir='C:\towerData\combined\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read Calibration Data from Site Details.xls or from mat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
radn_calib_xls = 'C:\towerData\SiteDetails.xlsx';
[n,t] = xlsread(radn_calib_xls,'CalibrationUpdates','A3:O100');
Calib_Headers = t(1,:);
Calib_Data_num = n(1:end,:);
Calib_Data_txt = t(2:end,:);

%Site info
SENSOR_SITE = Calib_Data_txt(:,1);
DATE_INSTALLED = Calib_Data_txt(:,5);

%calibration coefficients
CFACTOR1 = Calib_Data_num(:,5);
CFACTOR2 = Calib_Data_num(:,6);
CFACTOR3 = Calib_Data_num(:,7);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop through each site site level 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSite = sites2Proc
    %% LOAD FLUX DATA
    siteName = char(sites(iSite));
    inputfile = [mergedRootDir siteName '_MRG.mat'];

    disp(['loading :' inputfile]);
    eval(['load ' inputfile]);
    disp([inputfile ' loaded!']);
    
    HEADER = char(HMERGE);
    D = DMERGE;

    clear DMERGE;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set up filename for combined data
    combinedoutputfile = [combineRootDir siteName '_Combined.mat'];
    combinedoutputfilebak = [combineRootDir siteName '_Combined_' date '.mat'];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% What site are we working on
    disp(['Combining ' siteName]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%       Section 1: non-site specific adjustments and filtering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %NOTES
    %%%%%%These diagnostics check if irga and sonic are bad%%%%%%%%%

    % Red Alerts IRGA and SONIC out, but other sensors ok
    % Yellow Alerts, IRGA no good, SONIC and sensors ok

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DATA LOGGER DATA LOGGER DATA LOGGER DATA LOGGER DATA LOGGER DATA LOGGER
    %(1) Alerts from data logger
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %sonic and/or irga not working
    bad = D(143,:) < -998; %'Red_alert_Tot(1)'
    bad = bad | D(144,:) < -998; %'Yellow_alert_Tot(1)'
    bad = bad | D(146,:) < 3600; % 'IRGA_on_Tot' irga is on for half the possible measurements
    bad = bad | D(135,:) > 3600; %'irga_warning_Tot(1)' %first cut if 1/2 have warning on then remove
    bad = bad | D(135,:) == 480; %'irga_warning_Tot(1)' Hard 480s - no sure but maybe calibration ?
    bad = bad | D(135,:) > 3600; %'csat_warning_Tot(1)' out of 7200
    bad = bad | D(112,:) < minFlow(iSite);  %'irgaflow_Avg', l/min ==> need sufficient flow to reduce tube smearing
    bad = bad | D(110,:) < 50;  %'irgapress_Avg(1)'
    bad = bad | D(110,:) > 120; %pressure should be in this range 
        
    %NaN out covariances using sonic and irga
    D(63,bad)= NaN; %'Uz_co2_2(1)'
    D(64,bad)= NaN; %'Uz_h2o_2(1)'
    D(68,bad)= NaN; %'Ux_co2_2(1)'
    D(69,bad)= NaN; %'Ux_h2o_2(1)'
    D(72,bad)= NaN; %'Uy_co2_2(1)'
    D(73,bad)= NaN; %'Uy_h2o_2(1)'
    D(77,bad)= NaN; %'co2_Ts_2(1)' 
    D(79,bad)= NaN; %'h2o_Ts_2(1)'
    
    %Use this later to remove bad FAST data
    fast_bad_irga_sonic=bad;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %bad irga
    %IRGA not working
    bad = D(143,:) < -998; %'Red_alert_Tot(1)'
    bad = bad | D(144,:) < -998; %'Yellow_alert_Tot(1)'
    bad = bad | D(146,:) < 3600; % 'IRGA_on_Tot' irga is on for half the possible measurements
    bad = bad | D(135,:) > 3600; %'irga_warning_Tot(1)' %first cut if 1/2 have warning on then remove
    bad = bad | D(135,:) == 480; %'irga_warning_Tot(1)' Hard 480s - no sure but maybe calibration ?
    bad = bad | D(112,:) < minFlow(iSite);  %'irgaflow_Avg', l/min ==> need sufficient flow to reduce tube smearing
    bad = bad | D(110,:) < 50;  %'irgapress_Avg(1)'
    bad = bad | D(110,:) > 120; %pressure should be in this range 
    
    %NaN out irga
    D(75,bad)= NaN; %'co2_co2_2(1)'
    D(76,bad)= NaN; %'co2_h2o_2(1)'
    D(78,bad)= NaN; %'h2o_h2o_2(1)'
    D(85,bad)= NaN; %'co2_molar_Avg'
    D(86,bad)= NaN; %'h2o_molar_Avg'
    
    %Use this later to remove bad FAST data
    fast_bad_irga = bad;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SONIC is not working ==> Sonic fluxes are bad
    bad = D(143,:) < -998; %'Red_alert_Tot(1)'
    bad = bad | D(134,:) > 3600; %'csat_warning_Tot(1)' out of 7200
    
    %should NaN out all SONIC based measurements
    D(60:74,bad)= NaN; %'Uz_Uz_2(1) through 'Uy_Ts_2(1)'
    D(77,bad)= NaN; %'co2_Ts_2(1)'
    D(79:84,bad)= NaN; %'h2o_Ts_2(1) through'Ts_1_Avg(1)'

    %Use this later to remove bad FAST data
    fast_bad_sonic=bad;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %GOES GOES GOES GOES GOES GOES GOES GOES GOES GOES GOES GOES GOES GOES GOES
    %(2)Alerts from goes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %sonic and/or irga not working
    bad = D(269,:) < -998; %'GOES_red_alert'
    bad = bad | D(270,:) < -998; %'GOES_yellow_alert'
    bad = bad | D(272,:) < 3600; %GOES_irga_on
    bad = bad | D(261,:) > 3600; %'GOES_irga_warnings' %first cut if 1/2 have warning on then remove
    bad = bad | D(261,:) == 480; % calibration ?
    bad = bad | D(260,:) > 3600; %'GOES_csat_warnings'
    bad = bad | D(238,:) < minFlow(iSite); %'GOES_irgaflow_avg' l/min ==> need sufficient flow to reduce tube smearing
    bad = bad | D(236,:) < 50; %'GOES_irgapress_avg'
    bad = bad | D(236,:) > 120; %'GOES_irgapress_avg'

    %NaN out covariances using sonic and irga
    D(189,bad)= NaN; %'GOES_Uz_co2'
    D(190,bad)= NaN; %'GOES_Uz_h2o'
    D(194,bad)= NaN; %'GOES_Ux_co2'
    D(195,bad)= NaN; %'GOES_Ux_h2o'
    D(198,bad)= NaN; %'GOES_Uy_co2'
    D(199,bad)= NaN; %'GOES_Uy_h2o'
    D(203,bad)= NaN; %'GOES_co2_Ts' 
    D(205,bad)= NaN; %'GOES_h2o_Ts'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %bad Irga
    bad = D(269,:) < -998; %'GOES_red_alert'
    bad = bad | D(270,:) < -998; %'GOES_yellow_alert'
    bad = bad | D(272,:) < 3600; %GOES_irga_on
    bad = bad | D(261,:) > 3600; %'GOES_irga_warnings' %first cut if 1/2 have warning on then remove
    bad = bad | D(261,:) == 480; % calibration ?
    bad = bad | D(238,:) < minFlow(iSite); %'GOES_irgaflow_avg' l/min ==> need sufficient flow to reduce tube smearing
    bad = bad | D(236,:) < 50; %'GOES_irgapress_avg'
    bad = bad | D(236,:) > 120; %'GOES_irgapress_avg'

    %NaN out irga measurements
    D(201, bad)= NaN; %'GOES_co2_co2'
    D(202, bad)= NaN; %'GOES_co2_h2o'
    D(204, bad)= NaN; %'GOES_h2o_h2o'
    D(211, bad)= NaN; %'GOES_co2_avg'
    D(212, bad)= NaN; %'GOES_h2o_avg'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Now remove just bad sonic measurements
    bad = D(269,:) < -998; %'GOES_red_alert'
    bad = bad | D(260,:) > 3600; %'GOES_csat_warnings'
    
    %NaN out all Sonic based measurements
    D(186:200, bad)= NaN;
    D(203, bad)= NaN;
    D(205:210, bad)= NaN;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FAST %FAST  %FAST  %FAST  %FAST  %FAST  %FAST  %FAST  %FAST  %FAST 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %(3)FAST processed data
    %Assumption Alerts at data logger = alerts on fast processed data
    %==>NOTE: it is critical that the time stamp is being read correctly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %bad irga and sonic
    %D(29:30, fast_bad_irga_sonic)= NaN;
    %D(46:49, fast_bad_irga_sonic)= NaN;
    
    %bad irga
    %D(32:45, fast_bad_irga)= NaN;

    %bad SONIC measurements
    %D(5:28, fast_bad_sonic)= NaN;
    %D(31, fast_bad_sonic)= NaN;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Clean-up the Veg Index
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %data logger
    bad = D(92,:) <= 0; %VI sensors
    bad = bad | D(93,:) <= 0; %VI sensors
    bad = bad | D(94,:) <= 0; %VI sensors
    bad = bad | D(95,:) <= 0; %VI sensors
    bad = bad | D(96,:) <= 0; %VI sensors

    D(92:96, bad) = NaN; %bad NDVI
   
    %data logger
    bad = D(218,:) <= 0; %VI sensors
    bad = bad | D(219,:) <= 0; %VI sensors
    bad = bad | D(220,:) <= 0; %VI sensors
    bad = bad | D(221,:) <= 0; %VI sensors
    bad = bad | D(222,:) <= 0; %VI sensors
    D(218:222, bad) = NaN; %bad NDVI  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %the u_goes windspeed has to be corrected so we can do general calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The scaling is site dependant: Some windspeeds are ./100...others are .*100 
    %Scaling is for different time periods in different datasets. - awf

    %gain problem in the Ux direction
    if iSite==4; % LRgrass
        focus= D(1,:) > 529.252 & D(1,:) < 538.92;
        focus= focus | D(1,:) > 1256.938;
        D(207,~focus)= D(207,~focus).*100; %this is divided  
    elseif iSite == 5; %LRSage
        focus= D(1,:) > 529.22 & D(1,:) < 538.89;
        focus= focus | D(1,:) > 1243;
        D(207,~focus)= D(207,~focus).*100; %this is multiplied
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Most general ALL ALL ALL - entire dataset
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %bad over any measurement
    bad=D==-99999;
    bad = bad | D==-9999;
    bad = bad | D==-999;
    bad = bad | isempty(D) ==1;
    bad = bad | isfinite(D) == 0;
    D(bad)=NaN;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %***NOTES on general filtering***
    %removed all general filters - awf
        %there were 2 types
            %(1) multifilter - removed data using std, max, min, -999, etc.
            %(2) norm - homogenized different datasets before merging using
            %linear regression coefficents
        
    %***past comments***
    %(1)Have found problems with datalogger calculated diagnostics so
    % on 9/28/04 use only the goes data for the filtering

    %note 9/28/04   since we don't have good datalogger diagnostics at this
    %time use the goes data to filter the data logger data but when the
    %datalogger diagnostics are properly included into the merged files
    % then remove the code below whic assigns the datalogger iok flags to
    % the goes flags.

    %Data logger diagnositics typically agree with goes, although there are
    %a few departures. I used datalogger diagnostics independant of goes diagnostics. - awf
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Section 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Sensor Calibration --> turn to physical units
    
    %brought over code from Anne Kelly's calibrations - awf
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    where=char(sites(iSite));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %only Rn has been replaced 
    %%%add new changes manually
    DATE_INSTALLED = Calib_Data_txt(:,5);
    max_D=max(D(1,:))+1; % last day plus 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PAR_UP=PAR_In
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get gain
    Pin=[where '-PAR_UP'];
    good = strmatch(lower(Pin),lower(SENSOR_SITE));
    curr_multiplier = -1*CFACTOR1(good)*1000/604; % convert to current in amps 
    %V measured over a 604 ohm resistor (use V=IR) - checked with GW 4/6/12 - awf
    D(90,:)= D(90,:).* curr_multiplier; %dl 'par_in_Avg(1)'
    D(216,:)= D(216,:).* curr_multiplier; %'GOES_par_in_avg'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PAR_DOWN=PAR_Out
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get gain
    Pot=[where '-PAR_DN'];
    good = strmatch(lower(Pot),lower(SENSOR_SITE));
    curr_multiplier = -1*CFACTOR1(good)*1000/604; %convert to current in amps 604 Ohm resistor; 1000 accounts for mV ==> Amps
    D(91,:) = D(91,:) .* curr_multiplier; %dl
    D(217,:) = D(217,:) .* curr_multiplier; %'GOES_par_out_avg'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %pyrr up=pyrr in 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get gain
    Pyin=[where '-PYR_UP'];
    good = strmatch(lower(Pyin),lower(SENSOR_SITE));
    curr_multiplier = 1e-3*1/CFACTOR1(good); 
    D(88,:)= D(88,:).* curr_multiplier; %dl
    D(214,:)= D(214,:) .* curr_multiplier; %'GOES_pyrr_in_avg'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %pyrr down=pyrr out
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get gain
    Pyot=[where '-PYR_DN'];
    good = strmatch(lower(Pyot),lower(SENSOR_SITE));
    curr_multiplier = 1e-3*1/CFACTOR1(good); % we measure in mV - calibration is for 10^-6 V ==> 1e-3 (Check Manual)
    D(89,:)= D(89,:).* curr_multiplier; %dl
    D(215,:)= D(215,:) .* curr_multiplier; %'GOES_pyrr_in_avg'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %net radiation_pos
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %get gain 
    Rpos=[where '-RNET_POS'];
    good = strmatch(lower(Rpos),lower(SENSOR_SITE)); 
    c = CFACTOR1(good); 
    
    %- date subset
    curr_multiplier=D(1,:)*NaN; 
    DateChangeR=datenum(DATE_INSTALLED(good,1))-towerYearStart(iSite);
    good = D(1,:) >= 0 & D(1,:) < DateChangeR(2,1);
    curr_multiplier(good)=c(1);
    good = D(1,:) >= DateChangeR(2,1) & D(1,:) <  max_D;
    curr_multiplier(good)=c(2);

    %correct and replace data logger
    good = D(87,:) >=0;
    D(87,good)= D(87,good).* curr_multiplier(1,good);

    %correct and replace goes
    good = D(213,:) >=0;
    D(213,good)= D(213, good) .* curr_multiplier(1,good);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %net radiation_neg - sensor change is the same as net radiation_pos ==>
    %carry the same day through these calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get gain
    Rneg=[where '-RNET_NEG'];
    good = strmatch(lower(Rneg),lower(SENSOR_SITE));
    c = CFACTOR1(good); 
    
    %- date subset
    curr_multiplier=D(1,:)*NaN; 
    good = D(1,:) >= 0 & D(1,:) < DateChangeR(2,1);
    curr_multiplier(good)=c(1);
    good = D(1,:) >= DateChangeR(2,1) & D(1,:) <  max_D;
    curr_multiplier(good)=c(2);

    %correct and replace data logger
    good = D(87,:) < 0;
    D(87,good)= D(87,good).* curr_multiplier(1,good);

    %correct and replace goes
    good = D(213,:) < 0;
    D(213,good)= D(213, good) .* curr_multiplier(1,good);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Section 3:Calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %calculation of 27 new variables 
    %27 new variables will be added to the end of D

    %----------------------------------------------------------------
    %% Ustar
    %----------------------------------------------------------------
    %DATA LOGGER
    DL_Ustar = ( D(61,:).^2 + D(62,:).^2 ) .^ (1/4); %dl

    %Add to D and Header
    D(322,:)=DL_Ustar;
    HMERGE(322,1)={'Ustar_dl'};
    
    %goes
    GOES_Ustar   = ( D(187,:).^2 + D(188,:).^2 ) .^ (1/4);

    %Add to D and Header
    D(323,:)=GOES_Ustar;
    HMERGE(323,1)={'GOES_USTAR'};
    
    %----------------------------------------------------------------
    %% Wind Direction
    %----------------------------------------------------------------
     %The sonic orientation is measured from the N compase direction in 
    %a clockwise direction spoke with G Winston - awf
    
    %IMPORTANT NOTES:
    % all the processed data assumes that the sonics are pointing due east
    % (270°)
    
    % Actually the orientation varies from site to site  
    % Adjust for differences between sonic orientation
    %----------------------------------------------------------------
    %horizontal wind vectors
    %FAST
    u_mat=D(5,:);%'U mean, East to West'
    v_mat=D(6,:);%'V mean, North to South');
    
    %DATA LOGGER
    u_dl = D(81,:); %'Ux_1_Avg(1)');
    v_dl = D(82,:); %'Uy_1_Avg(1)');
    
    %goes
    v_goes = D(208,:); %'GOES_Uy_avg';
    u_goes = D(207,:); %'GOES_U_avg'; 
    %----------------------------------------------------------------
    %wind direction
    %FAST
    theta_mat = D(22,:);%'Wind Dir (met)');
    
    theta_mat = theta_mat - 90 + SonicOrientation(iSite);
    theta_mat(theta_mat>360)=theta_mat(theta_mat>360)- 360;
    theta_mat(theta_mat < 0)= theta_mat(theta_mat < 0)+360;

    %Add to D and Header
    D(324,:)=theta_mat;
    HMERGE(324,1)={'theta_mat'};
    
    %DATA LOGGER
    theta_dl   = (-atan2(v_dl,u_dl) + pi/2) * 180/pi;
    theta_dl(theta_dl < 0) = theta_dl(theta_dl < 0) + 360;
    
    %- 90 is because sonic data assumes that the sonics are pointed due east
    theta_dl   = theta_dl - 90 + SonicOrientation(iSite); 
    theta_dl(theta_dl>360)=theta_dl(theta_dl>360)- 360;
    theta_dl(theta_dl < 0)= theta_dl(theta_dl < 0)+360;

    %Add to D and Header
    D(325,:)=theta_dl;
    HMERGE(325,1)={'theta_dl'};
    
    %GOES
    theta_goes = (-atan2(v_goes,u_goes) + pi/2) *180/pi;
    theta_goes(theta_goes < 0)= theta_goes(theta_goes < 0)+360;
    
    theta_goes = theta_goes - 90 + SonicOrientation(iSite) ;
    theta_goes(theta_goes>360)=theta_goes(theta_goes>360)- 360;
    theta_goes(theta_goes < 0)= theta_goes(theta_goes < 0)+360;

    %Add to D and Header
    D(326,:)=theta_goes;
    HMERGE(326,1)={'theta_goes'};
    
    
    %-----------------------------------------------------------------
    %% Wind Speed
    %-----------------------------------------------------------------
    %FAST
    wspd_mat = sqrt(u_mat.^2 + v_mat.^2);
   
    %Add to D and Header
    D(327,:)=wspd_mat;
    HMERGE(327,1)={'wspd_mat'};
    
    %data logger 
    wspd_dl   = sqrt(u_dl.^2 + v_dl.^2);
    
    %Add to D and Header
    D(328,:)=wspd_dl;
    HMERGE(328,1)={'wspd_dl'};
    
    %GOES
    wspd_goes = sqrt(u_goes.^2 + v_goes.^2);
    
    %Add to D and Header
    D(329,:)=wspd_goes;
    HMERGE(329,1)={'wspd_goes'};
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      AIR DENSITY AND TRUE TEMP CALCULATIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Possible combinations for density
    %
    %   A: dry air density and dried sonic temperature RHO_DRY
    %   B: moist air density and dried sonic temperature
    %   C: dry air density and sonic temperature
    %   D: moist air density and sonic temperature
    %
    %     Use the Following Equations
    %         rhom_moist_air = p/RT
    %         PW = chi_h2o * rho_moist_air * RT
    %         Pa = p - PW
    %         rhom_dry_air = Pa /RT

    %----------------------------------------------------------------
    %%  True Temperature calculated by the Gaynor Equation
    %----------------------------------------------------------------
    % from another look an Sonic Anemometrty by Kaimal and Gaynor (1991) BLM
    % T = Ts / (1 + 0.32 e/P ), where T and Ts are in Kelvin and e and P
    % are in kPa
    %
    % equivalent to T = Ts / (1 - 0.00032 q ) where q is mmol of H2O per
    % mol air  
    
    %(This should be  + 0.00032q ???? I changed this.  -awf)

    %% (1) could fill here first to use optimal data for calculation
    %% (2) H2o concentrations are from irga - which are not calibrated well
    
    %GET SOME DATA FOR THE CALCULATION
    %sonic T
    Ts_dl   = D(84,:);% 'Ts_1_Avg(1)';
    Ts_goes = D(210 ,:);%'GOES_Ts_avg');
    Ts_mat = D(8,:);%'Ts');
    
    %!!!!!!!!!!!!!!!!CHECK THE UNITS!!!!!!!!!!!!!!!!!
    %IRGA water vapor concentrations
    h2o_dl   = D(86 ,:); %'h2o_molar_Avg');
    h2o_goes = D(212 ,:); %'GOES_h2o_avg');
    h2o_mat  = D(42 ,:); %'H2O Closed Path MEAN (mmol/mol)');
    
    %FAST
    Tkact_mat   = (Ts_mat + Tc)  ./ (1 + 0.00032 * h2o_mat );
    
    %Add to D and Header
    D(330,:)=Tkact_mat;
    HMERGE(330,1)={'Tkact_mat'};
    
    
    %DATALOGGER
    Tkact_dl    = (Ts_dl + Tc)   ./ (1 + 0.00032 * h2o_dl );
    
    %Add to D and Header
    D(331,:)=Tkact_dl;
    HMERGE(331,1)={'Tkact_dl'};
    
    %GOES
    Tkact_goes  = (Ts_goes + Tc) ./ (1 + 0.00032 * h2o_goes);
    
    %Add to D and Header
    D(332,:)=Tkact_goes;
    HMERGE(332,1)={'Tkact_goes'};
    
    
    %----------------------------------------------------------------
    %% Moist Air Density using dried Sonic Temperature
    %----------------------------------------------------------------
    %Press was oneAtmc;
    Press=sitePress(iSite); %kPa - awf
    
    %Density of Dried Air Reuqires Both the Sonic and Irga Data (molar
    %concentration - awf)
    rhom_ma_mat  = Press ./ (R_mol * Tkact_mat)  * 1e3; %units mol moist air per m3
    rhom_ma_dl   = Press ./ (R_mol * Tkact_dl)   * 1e3; %units mol moist air per m3
    rhom_ma_goes = Press ./ (R_mol * Tkact_goes) * 1e3; %units mol moist air per m3
    
    %Add to D and Header
    D(333,:)=rhom_ma_mat;
    HMERGE(333,1)={'rhom_ma_mat'};
    
    %Add to D and Header
    D(334,:)=rhom_ma_dl;
    HMERGE(334,1)={'rhom_ma_dl'};
    
    %Add to D and Header
    D(335,:)=rhom_ma_goes;
    HMERGE(335,1)={'rhom_ma_goes'};
    

    %----------------------------------------------------------------
    %% Partial pressure of water vapor
    %----------------------------------------------------------------
    Pw_mat  = rhom_ma_mat  .* h2o_mat  .* Tkact_mat  * R_mol * 1e-6; %kPa
    Pw_dl   = rhom_ma_dl   .* h2o_dl   .* Tkact_dl   * R_mol * 1e-6; %kPa
    Pw_goes = rhom_ma_goes .* h2o_goes .* Tkact_goes * R_mol * 1e-6; %kPa

    %----------------------------------------------------------------
    %% Partial pressure of Dry Air
    %---------------------------------------------------------------- 
    Pa_mat  = Press - Pw_mat; %kPa
    Pa_dl   = Press - Pw_dl;
    Pa_goes = Press - Pw_goes;

    %----------------------------------------------------------------
    %% Density of Dry Air
    %----------------------------------------------------------------
    rhom_da_mat  = 1e3 * Pa_mat  ./ (R_mol *  Tkact_mat); %mol/m^3
    rhom_da_dl   = 1e3 * Pa_dl   ./ (R_mol *  Tkact_dl);
    rhom_da_goes = 1e3 * Pa_goes ./ (R_mol *  Tkact_goes);

    %Add to D and Header
    D(336,:)=rhom_da_mat;
    HMERGE(336,1)={'rhom_da_mat'};
    
    %Add to D and Header
    D(337,:)=rhom_da_dl;
    HMERGE(337,1)={'rhom_da_dl'};
    
    %Add to D and Header
    D(338,:)=rhom_da_goes;
    HMERGE(338,1)={'rhom_da_goes'};
   
    %----------------------------------------------------------------
    %% (Water vapour molar density? - awf) Wet air molar density (mol/m3) and total mass density (kg/m3)
    %----------------------------------------------------------------
    rhom_water_mat = rhom_ma_mat - rhom_da_mat; %(mol/m3)
    rhom_water_dl   = rhom_ma_dl - rhom_da_dl;
    rhom_water_goes = rhom_ma_goes - rhom_da_goes;

    %total mass density (kg/m3)
    rhotot_mat = rhom_da_mat * Mw_da/1000  + rhom_water_mat * Mw_water/1000 ; %kg/m^3
    rhotot_dl   = rhom_da_dl * Mw_da/1000  + rhom_water_dl * Mw_water/1000 ;
    rhotot_goes = rhom_da_goes * Mw_da/1000  + rhom_water_goes * Mw_water/1000 ;
   
    %----------------------------------------------------------------
    %% Sensible Heat
    %----------------------------------------------------------------
    % get sensible heat flux from goes data
    % use the equation:    H = m_a * rho_a * CpAirc * <w'.T'>
    % where m_a = molar mass of air 28.966/1000
    % rho_a = density of air (mol m-3) calcd by rho = P/RT assum P is const (101.3 kPa) and R = 8.314 J/K/mol
    % CpAirc is the heat content of dry air (1004.67 J kg-1 K-1)
    % <w'.T'> is the 30 minute average of the vertical wind and the sonic
    % temperature
    
    %----------------------------------------------------------------
    %The correct way is to do this in the fastflux processing
    %estimate back calculation of SH
    %using recorded data
    sh_mat = D(29,:); %fast sensible heat
    
    Xw = D(42,:); %'H2O Closed Path MEAN (mmol/mol)'
    Cpv = CpH2Oc;

    %taken from the fastflux code
    Cpold = nanmean(CpAirc*(ones(size(Xw))-Xw) + Cpv.*Xw);
    Cpnew = nanmean(CpAirc*(ones(size(Xw))-(1e-3*Xw)) + Cpv.*(1e-3*Xw));

    mult=Cpnew/Cpold;
    
    D(29,:)=mult .* sh_mat; %adjust for new Cp
    
    % declare some constants initially
    m_a= Mw_da / 1000; %units of kg dry air / mol
    %----------------------------------------------------------------
    
    %DATA LOGGER
    sh_dl   = m_a * CpAirc * rhom_da_mat.* D(65,:); %units of W m-2; %D(65.:)='Uz_Ts_2(1)';

    %Add to D and Header
    D(339,:)= sh_dl;
    HMERGE(339,1)={'sh_dl W m-2'};
    
    
    %GOES
    sh_goes = m_a * CpAirc * rhom_da_goes .* D(191,:); %units of W m-2;

    %Add to D and Header
    D(340,:)= sh_goes;
    HMERGE(340,1)={'sh_goes W m-2'};
    
    %----------------------------------------------------------------
    %% Approximate Sensible Heat based purely on Sonic Signals %%%%%
    %----------------------------------------------------------------
    h_buoyant_dl   = Mw_da/1000 * 38.6 * CpAirc * D(65,:);
    
    %Add to D and Header
    D(341,:)= h_buoyant_dl ;
    HMERGE(341,1)={'h_buoyant_dl '};
    
    h_buoyant_goes = Mw_da/1000 * 38.6 * CpAirc * D(191,:);

    %Add to D and Header
    D(342,:)= h_buoyant_goes;
    HMERGE(342,1)={'h_buoyant_goes'};
    
    %----------------------------------------------------------------
    %% Latent Heat
    %----------------------------------------------------------------

    % %for the slow data use the equation
    % %
    % %    LH = m_w / 1000 * rho_a * Lv * <w'q'> /1000 :units=W/m2
    % %    m_w = molecular weight of water = 18 g/mol
    % %    rho_a = density of air (mol/m3)
    % %    Lv = Latent heat of vaporization (J/kg)
    % %    <w'q'> =covaraiance of w (m/s) and q (mmol/mol)
    % %
    
    %LOGGER
    %calculate Latent heat of vaporization using the sonic temperature 
    Lv_dl   = (2.501-0.00237.*(Ts_dl-Tc))*10^6; %units: J/kg
    
    %calculate latent heat
    lh_dl   = 1e-6 * 18 * rhom_da_dl.*Lv_dl .* D(64,:); %'Uz_h2o_2(1)'
    
    %Add to D and Header
    D(343,:)=lh_dl;
    HMERGE(343,1)={'lh_dl'};
    
    %GOES
    %calculate Latent heat of vaporization using the sonic temperature
    Lv_goes = (2.501-0.00237.*(Ts_goes-Tc))*10^6; %units: J/kg
    
    %calculate latent heat
    lh_goes = 1e-6 * 18 * rhom_da_goes.*Lv_goes .* D(190,:); %'GOES_Uz_h2o'

    %Add to D and Header
    D(344,:)=lh_goes;
    HMERGE(344,1)={'lh_goes'};
    
    %pressure correction????
    
    %could use shaw model gains to correct for tube smearing here
    %     %now correct for the delay and smearing in the datalogger and goes data using the Shaw Tau model
    %     gain_fh2o = chset( D,HEADER,'Gain for h2o flux using sonic temperature flux and Shaw Model tau');
    %     %but only use gain if it is a factor between 0.5 and 1.5;
    %
    %     lh_dl   = lh_dl.*gain_fh2o;
    %     lh_goes = lh_goes.*gain_fh2o;
    
    %----------------------------------------------------------------
    %% CO2 fluxes
    %----------------------------------------------------------------
    %UNITS: mol dry air /m^3 * m/s * micromol CO2/ mol air ==> micromol CO2 / m^2 / s
    
    %datalogger
    Fco2_dl = rhom_da_dl.* D(63,:);%'Uz_co2_2(1)'; 
    
    %save the data
    D(345,:)=Fco2_dl;
    HMERGE(345,1)={'Fco2_dl'};
    
    %goes
    Fco2_goes  = rhom_da_goes.* D(189,:); %'GOES_Uz_co2';
    
    %save the data
    D(346,:)=Fco2_goes;
    HMERGE(346,1)={'Fco2_goes'};
    
    %----------------------------------------------------------------
    %% H2O fluxes
    %----------------------------------------------------------------
    %data logger
    Fh2o_dl = rhom_da_dl.* D(64,:);%'Uz_h2o_2(1)' 
    
    %save the data
    D(347,:)=Fh2o_dl;
    HMERGE(347,1)={'Fh2o_dl'};
    
    %goes
    Fh2o_goes = rhom_da_goes.* D(190,:);%'GOES_Uz_h2o'
    
    %save the data
    D(348,:)=Fh2o_goes;
    HMERGE(348,1)={'Fh2o_goes'};
    %now correct for the delay and smearing in the datalogger and goes data using the Shaw Tau model

    %     fh2o_dl   = fh2o_dl.*gain_fh2o;
    %     fh2o_goes = fh2o_goes.*gain_fh2o;
    
    %----------------------------------------------------------------------
    %%   Water Vapor Pressure (HMP Probe)
    %----------------------------------------------------------------------
    %Ideal Gas Law
    %backup density when irga is not working
    rhom_ma_wetsonic_dl   = Press ./ (R_mol * (Ts_dl   + Tc)) * 1e3; %units mol moist air per m3
    rhom_ma_wetsonic_goes = Press ./ (R_mol * (Ts_goes + Tc)) * 1e3; %units mol moist air per m3

    %DATA LOGGER
    %D(100,:)= 'h2o_hmp_Avg' are in units of g/m3 (h2o)
    %convert to units millimol h2o/mol moist air
    h2o_hmp_dl=(D(100,:)./ rhom_ma_dl) * (1/18) * 1000 ;
    D(100,:)=h2o_hmp_dl;
    
    %Think about this
    %h2o_hmp_dl_noIRGA=(h2o_hmp_dl./ rhom_ma_wetsonic_dl) * (1/18) * 1000 ;
    
    %GOES
    %convert to units millimol h2o/mol moist air
    h2o_hmp_goes=(D(226,:)./ rhom_ma_goes) * (1/18) * 1000 ;
    D(226,:)= h2o_hmp_goes;
    
    %Think about this
    %h2o_hmp_goes_noIRGA=(h2o_hmp_goes./ rhom_ma_wetsonic_goes) * (1/18) * 1000 ;
    
    %----------------------------------------------------------------------
    %redo HEADER because we added new rows to the data D
    HEADER = char(HMERGE);
    %----------------------------------------------------------------------
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Section 4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Site specific filtering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
        if iSite == 1
            %go to this subdirectory
        elseif iSite ==2
            %go to this subdirectory
        elseif iSite ==3
            %go to this subdirectory
        elseif iSite ==4 %LR grassland
            D=Site_specific_LR_Grass2(HEADER, D);
        elseif iSite ==5
            %go to this subdirectory
            D=Site_specific_LR_Sage3(HEADER, D);
        elseif iSite ==6
            %go to this subdirectory
        elseif iSite ==7
            %go to this subdirectory
        elseif iSite ==8
            %go to this subdirectory
        elseif iSite ==9
            %go to this subdirectory
        else
            disp(['you have no site specific filters for ' siteName]) 
            disp('you may have introduced bad data into your output')
            pause
        end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Section 5    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Here is where we combined
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %Each dataset has been filtered and corrected for sensor calibrations 
    %Now we will combined these different data streams into 1 data table
    %We combined data streams in the following order: 1) fast, 2)datalogger, 3)goes 
    
    %What is the variable? orginized by column
    Cleaned_Header={};
    
    %Here is the data for the variable in Cleaned_Header
    ncols=size(D,2);
    Cleaned_D=ones(29, ncols)*NaN;
    
    %What data stream did the data come from?
    DataSource=ones(29, ncols)*NaN; %1 = processed; 2=datalogger; 3=GOES
    
    %----------------------------------------------------------------------
    %1.Time
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(1,:)= D(1,:);
    
    %Fil Header
    Cleaned_Header(1,1)={'TIME'};
    
    %----------------------------------------------------------------------
    %2.Tsonic
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(2,:)=D(8,:);%fast
    bad=isnan(Cleaned_D(2,:));
    Cleaned_D(2,bad)=D(84,bad);%dl
    bad=isnan(Cleaned_D(2,:));
    Cleaned_D(2,bad)=D(210,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(8,:));
    good_dl=(isnan(D(8,:)) & ~isnan(D(84,:)));
    good_goes=(isnan(D(8,:)) & isnan(D(84,:)) & ~isnan(D(210,:)));
    DataSource(2,good_mat)=1;
    DataSource(2,good_dl)=2;
    DataSource(2,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(2,1)={'TSONIC'};
    
    %----------------------------------------------------------------------
    %3.Tactual
    %----------------------------------------------------------------------
    %****************************************
    %Fill Data
    Cleaned_D(3,:)=D(330,:);%fast
    bad=isnan(Cleaned_D(3,:));
    Cleaned_D(3,bad)=D(331,bad);%dl
    bad=isnan(Cleaned_D(3,:));
    Cleaned_D(3,bad)=D(332,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(330,:));
    good_dl=(isnan(D(330,:)) & ~isnan(D(331,:)));
    good_goes=(isnan(D(330,:)) & isnan(D(331,:)) & ~isnan(D(332,:)));
    DataSource(3,good_mat)=1;
    DataSource(3,good_dl)=2;
    DataSource(3,good_goes)=3;
    
    
    %Fill Header
    Cleaned_Header(3,1)={'TACTUAL'};
    %****************************************
    %----------------------------------------------------------------------
    %4.Ustar
    %----------------------------------------------------------------------
    %****************************************
    %Fill Data
    Cleaned_D(4,:)=D(21,:);
    bad=isnan(Cleaned_D(4,:));
    Cleaned_D(4,bad)=D(322,bad);
    bad=isnan(Cleaned_D(4,:));
    Cleaned_D(4,bad)=D(323,bad);
    
    %Where is the data from?
    good_mat=~isnan(D(21,:));
    good_dl=(isnan(D(21,:)) & ~isnan(D(322,:)));
    good_goes=(isnan(D(21,:)) & isnan(D(322,:)) & ~isnan(D(323,:)));
    DataSource(4,good_mat)=1;
    DataSource(4,good_dl)=2;
    DataSource(4,good_goes)=3;
    
    %put USTAR source into Cleaned_D
    Cleaned_D(25,:)=DataSource(4,:);
    bad=Cleaned_D(25,:)==0;
    Cleaned_D(25,bad)=NaN;
  
    %Fill Header
    Cleaned_Header(4,1)={'USTAR'};
    Cleaned_Header(25,1)={'USTAR_SOURCE_FLAG'};
    %****************************************
    %----------------------------------------------------------------------
    %5.Wind Direction
    %----------------------------------------------------------------------
    %****************************************
    %Fill Data
    Cleaned_D(5,:)=D(324,:);%fast
    bad=isnan(Cleaned_D(5,:));
    Cleaned_D(5,bad)=D(325,bad);%dl
    bad=isnan(Cleaned_D(5,:));
    Cleaned_D(5,bad)=D(326,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(324,:));
    good_dl=(isnan(D(324,:)) & ~isnan(D(325,:)));
    good_goes=(isnan(D(324,:)) & isnan(D(325,:)) & ~isnan(D(326,:)));
    DataSource(5,good_mat)=1;
    DataSource(5,good_dl)=2;
    DataSource(5,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(5,1)={'WDIR'};
    %****************************************
    %----------------------------------------------------------------------
    %6.Wind Speed
    %----------------------------------------------------------------------
    %****************************************
    %Fill Data
    Cleaned_D(6,:)=D(327,:);%fast
    bad=isnan(Cleaned_D(6,:));
    Cleaned_D(6,bad)=D(328,bad);%dl
    bad=isnan(Cleaned_D(6,:));
    Cleaned_D(6,bad)=D(329,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(327,:));
    good_dl=(isnan(D(327,:)) & ~isnan(D(328,:)));
    good_goes=(isnan(D(327,:)) & isnan(D(328,:)) & ~isnan(D(329,:)));
    DataSource(6,good_mat)=1;
    DataSource(6,good_dl)=2;
    DataSource(6,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(6,1)={'WSPD'};
    %****************************************
    %----------------------------------------------------------------------
    %7. Dry Density
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(7,:)=D(336,:);%fast
    bad=isnan(Cleaned_D(7,:));
    Cleaned_D(7,bad)=D(337,bad);%dl
    bad=isnan(Cleaned_D(7,:));
    Cleaned_D(7,bad)=D(338,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(336,:));
    good_dl=(isnan(D(336,:)) & ~isnan(D(337,:)));
    good_goes=(isnan(D(336,:)) & isnan(D(337,:)) & ~isnan(D(338,:)));
    DataSource(7,good_mat)=1;
    DataSource(7,good_dl)=2;
    DataSource(7,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(7,1)={'DENSITY_DRY'};
    %****************************************
    %----------------------------------------------------------------------
    %8. Moist Density
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(8,:)=D(333,:);%fast
    bad=isnan(Cleaned_D(8,:));
    Cleaned_D(8,bad)=D(334,bad);%dl
    bad=isnan(Cleaned_D(8,:));
    Cleaned_D(8,bad)=D(335,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(333,:));
    good_dl=(isnan(D(333,:)) & ~isnan(D(334,:)));
    good_goes=(isnan(D(333,:)) & isnan(D(334,:)) & ~isnan(D(335,:)));
    DataSource(8,good_mat)=1;
    DataSource(8,good_dl)=2;
    DataSource(8,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(8,1)={'DENSITY_MOIST'};
    %----------------------------------------------------------------------
    %9. Sensible heat
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(9,:)=D(29,:);%fast
    bad=isnan(Cleaned_D(9,:));
    Cleaned_D(9,bad)=D(339,bad);%dl
    bad=isnan(Cleaned_D(9,:));
    Cleaned_D(9,bad)=D(340,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(29,:));
    good_dl=(isnan(D(29,:)) & ~isnan(D(339,:)));
    good_goes=(isnan(D(29,:)) & isnan(D(339,:)) & ~isnan(D(340,:)));
    DataSource(9,good_mat)=1;
    DataSource(9,good_dl)=2;
    DataSource(9,good_goes)=3;
    
    %put source into Cleaned_D
    Cleaned_D(26,:)=DataSource(9,:);
    bad=Cleaned_D(26,:)==0;
    Cleaned_D(26,bad)=NaN;
    
    %Fill Header
    Cleaned_Header(9,1)={'SENSIBLE_HEAT'};
    Cleaned_Header(26,1)={'SENSIBLE_HEAT_SOURCE_FLAG'};
    %----------------------------------------------------------------------
    %10. Latent heat
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(10,:)=D(30,:);%fast
    bad=isnan(Cleaned_D(10,:));
    Cleaned_D(10,bad)=D(343,bad);%dl
    bad=isnan(Cleaned_D(10,:));
    Cleaned_D(10,bad)=D(344,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(30,:));
    good_dl=(isnan(D(30,:)) & ~isnan(D(343,:)));
    good_goes=(isnan(D(30,:)) & isnan(D(343,:)) & ~isnan(D(344,:)));
    DataSource(10,good_mat)=1;
    DataSource(10,good_dl)=2;
    DataSource(10,good_goes)=3;
    
    %put source into Cleaned_D
    Cleaned_D(27,:)=DataSource(10,:);
    bad=Cleaned_D(27,:)==0;
    Cleaned_D(27,bad)=NaN;
    
    %Fill Header
    Cleaned_Header(10,1)={'LATENT_HEAT'};
    Cleaned_Header(27,1)={'LATENT_HEAT_SOURCE_FLAG'};
    %----------------------------------------------------------------------
    %11. CO2 concentration
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(11,:)=D(35,:);%fast
    bad=isnan(Cleaned_D(11,:));
    Cleaned_D(11,bad)=D(85,bad);%dl
    bad=isnan(Cleaned_D(11,:));
    Cleaned_D(11,bad)=D(211,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(35,:));
    good_dl=(isnan(D(35,:)) & ~isnan(D(85,:)));
    good_goes=(isnan(D(35,:)) & isnan(D(85,:)) & ~isnan(D(211,:)));
    DataSource(11,good_mat)=1;
    DataSource(11,good_dl)=2;
    DataSource(11,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(11,1)={'CO2_CONC'};
    %----------------------------------------------------------------------
    %12. H2O concentration
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(12,:)=D(42,:);%fast
    bad=isnan(Cleaned_D(12,:));
    Cleaned_D(12,bad)=D(86,bad);%dl
    bad=isnan(Cleaned_D(12,:));
    Cleaned_D(12,bad)=D(212,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(42,:));
    good_dl=(isnan(D(42,:)) & ~isnan(D(86,:)));
    good_goes=(isnan(D(42,:)) & isnan(D(86,:)) & ~isnan(D(212,:)));
    DataSource(12,good_mat)=1;
    DataSource(12,good_dl)=2;
    DataSource(12,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(12,1)={'H2O_CONC'};
    %----------------------------------------------------------------------
    %13. FCO2
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(13,:)=D(46,:);%fast
    bad=isnan(Cleaned_D(13,:));
    Cleaned_D(13,bad)=D(345,bad);%dl
    bad=isnan(Cleaned_D(13,:));
    Cleaned_D(13,bad)=D(346,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(46,:));
    good_dl=(isnan(D(46,:)) & ~isnan(D(345,:)));
    good_goes=(isnan(D(46,:)) & isnan(D(345,:)) & ~isnan(D(346,:)));
    DataSource(13,good_mat)=1;
    DataSource(13,good_dl)=2;
    DataSource(13,good_goes)=3;
    
    %put source into Cleaned_D
    Cleaned_D(28,:)=DataSource(13,:);
    bad=Cleaned_D(28,:)==0;
    Cleaned_D(28,bad)=NaN;
    
    %Fill Header
    Cleaned_Header(13,1)={'CO2_FLUX'};
    Cleaned_Header(28,1)={'CO2_FLUX_SOURCE_FLAG'};
    %----------------------------------------------------------------------
    %14. FH2O
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(14,:)=D(47,:);%fast
    bad=isnan(Cleaned_D(14,:));
    Cleaned_D(14,bad)=D(347,bad);%dl
    bad=isnan(Cleaned_D(14,:));
    Cleaned_D(14,bad)=D(348,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(47,:));
    good_dl=(isnan(D(47,:)) & ~isnan(D(347,:)));
    good_goes=(isnan(D(47,:)) & isnan(D(347,:)) & ~isnan(D(348,:)));
    DataSource(14,good_mat)=1;
    DataSource(14,good_dl)=2;
    DataSource(14,good_goes)=3;
    
    %put source into Cleaned_D
    Cleaned_D(29,:)=DataSource(14,:);
    bad=Cleaned_D(29,:)==0;
    Cleaned_D(29,bad)=NaN;
    
    %Fill Header
    Cleaned_Header(14,1)={'H2O_FLUX'};
    Cleaned_Header(29,1)={'H2O_FLUX_SOURCE_FLAG'};
    %----------------------------------------------------------------------
    %15. Rn
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(15,:)=D(87,:);%dl
    bad=isnan(Cleaned_D(15,:));
    Cleaned_D(15,bad)=D(213,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(87,:));
    good_goes=(isnan(D(87,:)) & ~isnan(D(213,:)));
    DataSource(15,good_dl)=2;
    DataSource(15,good_goes)=3;

    %Fill Header
    Cleaned_Header(15,1)={'RNET'};
    %----------------------------------------------------------------------
    %16. PAR_In
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(16,:)=D(90,:);%dl
    bad=isnan(Cleaned_D(16,:));
    Cleaned_D(16,bad)=D(216,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(90,:));
    good_goes=(isnan(D(90,:)) & ~isnan(D(216,:)));
    DataSource(16,good_dl)=2;
    DataSource(16,good_goes)=3;

    %Fill Header
    Cleaned_Header(16,1)={'PAR_IN'};
    %----------------------------------------------------------------------
    %17. PAR_OUT
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(17,:)=D(91,:);%dl
    bad=isnan(Cleaned_D(17,:));
    Cleaned_D(17,bad)=D(217,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(91,:));
    good_goes=(isnan(D(91,:)) & ~isnan(D(217,:)));
    DataSource(17,good_dl)=2;
    DataSource(17,good_goes)=3;

    %Fill Header
    Cleaned_Header(17,1)={'PAR_OUT'};

    %----------------------------------------------------------------------
    %18. Incoming solar radiation - pyranometer = SOLAR_IN
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(18,:)=D(88,:);%dl
    bad=isnan(Cleaned_D(18,:));
    Cleaned_D(18,bad)=D(214,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(88,:));
    good_goes=(isnan(D(88,:)) & ~isnan(D(214,:)));
    DataSource(18,good_dl)=2;
    DataSource(18,good_goes)=3;

    %Fill Header
    Cleaned_Header(18,1)={'SOLAR_IN'};

    %----------------------------------------------------------------------
    %19. Outgoing solar radiation - pyranometer = SOLAR_OUT
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(19,:)=D(89,:);%dl
    bad=isnan(Cleaned_D(19,:));
    Cleaned_D(19,bad)=D(215,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(89,:));
    good_goes=(isnan(D(89,:)) & ~isnan(D(215,:)));
    DataSource(19,good_dl)=2;
    DataSource(19,good_goes)=3;

    %Fill Header
    Cleaned_Header(19,1)={'SOLAR_OUT'};
    
    %----------------------------------------------------------------------
    %20. HMP_Temp 
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(20,:)=D(97,:);%dl
    bad=isnan(Cleaned_D(20,:));
    Cleaned_D(20,bad)=D(223,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(97,:));
    good_goes=(isnan(D(97,:)) & ~isnan(D(223,:)));
    DataSource(20,good_dl)=2;
    DataSource(20,good_goes)=3;

    %Fill Header
    Cleaned_Header(20,1)={'T_HMP'}; 
    
    %----------------------------------------------------------------------
    %21. HMP_RH
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(21,:)=D(98,:);%dl
    bad=isnan(Cleaned_D(21,:));
    Cleaned_D(21,bad)=D(224,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(98,:));
    good_goes=(isnan(D(98,:)) & ~isnan(D(224,:)));
    DataSource(21,good_dl)=2;
    DataSource(21,good_goes)=3;

    %Fill Header
    Cleaned_Header(21,1)={'RH'}; 
        
    %----------------------------------------------------------------------
    %22. T_107
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(22,:)=D(102,:);%dl
    bad=isnan(Cleaned_D(22,:));
    Cleaned_D(22,bad)=D(228,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(102,:));
    good_goes=(isnan(D(102,:)) & ~isnan(D(228,:)));
    DataSource(22,good_dl)=2;
    DataSource(22,good_goes)=3;

    %Fill Header
    Cleaned_Header(22,1)={'T_107'}; 
    
    %----------------------------------------------------------------------
    %23. Rain
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(23,:)=D(103,:);%dl
    bad=isnan(Cleaned_D(23,:));
    Cleaned_D(23,bad)=D(229,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(103,:));
    good_goes=(isnan(D(103,:)) & ~isnan(D(229,:)));
    DataSource(23,good_dl)=2;
    DataSource(23,good_goes)=3;

    %Fill Header
    Cleaned_Header(23,1)={'RAIN'}; 
    
    %----------------------------------------------------------------------
    %24. Water vapor concentration
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(24,:)=D(100,:);%dl
    bad=isnan(Cleaned_D(24,:));
    Cleaned_D(24,bad)=D(226,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(100,:));
    good_goes=(isnan(D(100,:)) & ~isnan(D(226,:)));
    DataSource(24,good_dl)=2;
    DataSource(24,good_goes)=3;

    %Fill Header
    Cleaned_Header(24,1)={'H2O_HMP'}; 
    
    %----------------------------------------------------------------------
    %30. NDVI
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(30,:)=D(96,:);%dl
    bad=isnan(Cleaned_D(30,:));
    Cleaned_D(30,bad)=D(222,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(96,:));
    good_goes=(isnan(D(96,:)) & ~isnan(D(222,:)));
    DataSource(30,good_dl)=2;
    DataSource(30,good_goes)=3;

    %Fill Header
    Cleaned_Header(30,1)={'NDVI'}; 
    
    %----------------------------------------------------------------------
    %31. RED IN
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(31,:)=D(92,:);%dl
    bad=isnan(Cleaned_D(31,:));
    Cleaned_D(31,bad)=D(218,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(92,:));
    good_goes=(isnan(D(92,:)) & ~isnan(D(218,:)));
    DataSource(31,good_dl)=2;
    DataSource(31,good_goes)=3;

    %Fill Header
    Cleaned_Header(31,1)={'Red In'}; 
    
    %----------------------------------------------------------------------
    %32. RED OUT
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(32,:)=D(94,:);%dl
    bad=isnan(Cleaned_D(32,:));
    Cleaned_D(32,bad)=D(220,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(94,:));
    good_goes=(isnan(D(94,:)) & ~isnan(D(220,:)));
    DataSource(32,good_dl)=2;
    DataSource(32,good_goes)=3;

    %Fill Header
    Cleaned_Header(32,1)={'Red Out'}; 
    
    %----------------------------------------------------------------------
    %33. NIR In
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(33,:)=D(93,:);%dl
    bad=isnan(Cleaned_D(33,:));
    Cleaned_D(33,bad)=D(219,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(93,:));
    good_goes=(isnan(D(93,:)) & ~isnan(D(219,:)));
    DataSource(33,good_dl)=2;
    DataSource(33,good_goes)=3;

    %Fill Header
    Cleaned_Header(33,1)={'NIR In'}; 
    
    %----------------------------------------------------------------------
    %34. NIR Out
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(34,:)=D(95,:);%dl
    bad=isnan(Cleaned_D(34,:));
    Cleaned_D(34,bad)=D(221,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(95,:));
    good_goes=(isnan(D(95,:)) & ~isnan(D(221,:)));
    DataSource(34,good_dl)=2;
    DataSource(34,good_goes)=3;

    %Fill Header
    Cleaned_Header(34,1)={'NIR Out'}; 
    
    %----------------------------------------------------------------------
    %35. Soil Temperature and Moisture 
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(35:44,:)=D(284:293,:);%dl
    Cleaned_D(45:71,:)=D(295:321,:);%dl
    
    %Where is the data from?
    DataSource(35:71,:)=NaN;

    %Fill Header
    Cleaned_Header(35,1)={'Soil_M(1)'}; 
    Cleaned_Header(36,1)={'Soil_M(2)'}; 
    Cleaned_Header(37,1)={'Soil_M(3)'}; 
    Cleaned_Header(38,1)={'Soil_M(4)'}; 
    Cleaned_Header(39,1)={'Fuel_M(1)'}; 
    Cleaned_Header(40,1)={'Fuel_M(2)'}; 
    Cleaned_Header(41,1)={'Soil_T(1)'}; 
    Cleaned_Header(42,1)={'Soil_T(2)'}; 
    Cleaned_Header(43,1)={'Soil_T(3)'}; 
    Cleaned_Header(44,1)={'Soil_T(4)'}; 
    Cleaned_Header(45,1)={'LWS(1)'}; 
    Cleaned_Header(46,1)={'LWS(2)'}; 
    Cleaned_Header(47,1)={'LWS(3)'}; 
    Cleaned_Header(48,1)={'StartT_C(1)'}; 
    Cleaned_Header(49,1)={'StartT_C(2)'}; 
    Cleaned_Header(50,1)={'StartT_C(3)'}; 
    Cleaned_Header(51,1)={'StartT_C(4)'}; 
    Cleaned_Header(52,1)={'StartT_C(5)'}; 
    Cleaned_Header(53,1)={'StartT_C(6)'}; 
    Cleaned_Header(54,1)={'T_1s_C(1)'}; 
    Cleaned_Header(55,1)={'T_1s_C(2)'}; 
    Cleaned_Header(56,1)={'T_1s_C(3)'}; 
    Cleaned_Header(57,1)={'T_1s_C(4)'}; 
    Cleaned_Header(58,1)={'T_1s_C(5)'}; 
    Cleaned_Header(59,1)={'T_1s_C(6)'}; 
    Cleaned_Header(60,1)={'T_30s_C(1)'}; 
    Cleaned_Header(61,1)={'T_30s_C(2)'}; 
    Cleaned_Header(62,1)={'T_30s_C(3)'}; 
    Cleaned_Header(63,1)={'T_30s_C(4)'}; 
    Cleaned_Header(64,1)={'T_30s_C(5)'}; 
    Cleaned_Header(65,1)={'T_30s_C(6)'}; 
    Cleaned_Header(66,1)={'DelT_C(1)'}; 
    Cleaned_Header(67,1)={'DelT_C(2)'}; 
    Cleaned_Header(68,1)={'DelT_C(3)'}; 
    Cleaned_Header(69,1)={'DelT_C(4)'}; 
    Cleaned_Header(70,1)={'DelT_C(5)'}; 
    Cleaned_Header(71,1)={'DelT_C(6)'}; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Section 6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Flip to maintain format
    Cleaned_D=Cleaned_D';
    DataSource=DataSource';
    
    %% Save
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

fnamecleaned = [combineRootDir  siteName '_cleaned'];

save(fnamecleaned, 'Cleaned_D', 'Cleaned_Header', 'DataSource');

end
