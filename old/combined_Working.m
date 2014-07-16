function combined_Working(sites2Proc, saveon, combinInfo)
%Inputs
%What sites do you want to process? see var_defs sites table to relate
%number to site

%saveon=1 saves the combined data
%combinInfo =1 runs and saves ancillary information on what you will
%combined

%Output 
% this output data is saved in 'C:towerData\cleaned' 
%'Cleaned_D', - combined best dataset 
%'DataSource', - what data stream 
%'Cleaned_Header', - header for Cleaned_D and DataSource

%'D', - all data - merge with some additional calculations
%'HMERGE' - D's header

%this output data is saved in 'C:towerData\cleaned'
%D_obs_out will give the observations in D that were filtered out 
%N_obs - summary of omiited observations
%headerN_obs - N_obs header
%information on date, high freq Gain, and delay
%--------------------------------------------------------------------------
%origin:
%Nearly all of this code comes from the Anne Kelly's combined_new program circa
%4/3/12
%Reorginization of code by awf %summer 2012
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
%There are 7 sections
%(1)non-site specific filtering
%(2)Sensor Calibration --> turn to physical units
%(3)Site specific filtering - uses additional site specific programs
%(4)Calculations
%(5)Mop-up - remove outliers from calculations
%(6)homogenize and combined data streams
%(7)save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Initialilizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Site-specific info
global sites sitePress towerYearStart iSite SonicOrientation minFlow minPress
% Paths etc.
path(path, 'C:\towerData\processingScripts\subroutines');
% flux processing parameters
global R_mol Tc CpAirc Mw_da Mw_water
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
    
    if combinInfo ==1
        N_obs=[];
        %D_obs_out will give the observations in D that were omitted 
        size_D=size(D, 2);
        D_obs_out=ones(355,size_D).*NaN; %this will get big
        
       [N_obs, headerN_obs, D_obs_out] = combinedinfo(siteName, D, D_obs_out, N_obs, 1);
    end
    %Take this out - it is rigged for prelim Sierra Site check
    %n=D(1,:)*NaN;
    %D2=[D(1:107,:); n; D(108:112,:); n; D(113:114,:); n; D(115:168,:); n; D(169:172,:); n;n;n;n; D(173:225,:); n; D(226:230,:); n; D(231:232,:); n; D(233:264,:)];
    %D2=[D(1:107,:); n; D(108:112,:); n; D(113:114,:); n; D(115:167,:); n;n;n;n;n;n;n;n;n;n;D(168:204,:); n;n;n;n;n; D(205:215,:); n; D(216:220,:); n; D(221:222,:); n; D(223:257,:)]; % old shorthair
    %D2(5:180,2:44546)=D2(5:180,1:44545); % old shorthair
    %D2(5:180,1)=NaN; % old shorthair
    %D=D2;
    %end take out 
   
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

    %use dl alerts to NaN fast processed data
    
    %take a conservative approach - if the alert is missing ==> omit
    %measurements
    
    %checked logic on alerts - looks ok -awf 10/4/2012
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DATA LOGGER DATA LOGGER DATA LOGGER DATA LOGGER DATA LOGGER DATA LOGGER
    %(1) Alerts from data logger
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %'csat_warning and irga warnings - use this to remove -999s - important
    %for filtering
    
    bad = D(134,:) < -1; 
    D(134,bad)=NaN;
    
    bad = D(135,:) < -1; 
    D(135,bad)=NaN;
    
    %Removed flow filter because the flow meter was missing periods or
    %giving zeros, which caused us to remove good data.  Flow and pressure are related at most sites
    %==> used the pressure filter instead of flow. - awf
    
    %sonic and/or irga not working ==> remove only sonic and irga based
    %measurements
    good = D(143,:) > -998; %'Red_alert_Tot(1)'
    good = good & D(144,:) > -998; %'Yellow_alert_Tot(1)'
    good = good & D(146,:) > 3600; % 'IRGA_on_Tot' irga is on for half the possible measurements
    good = good & D(135,:) < 3600; %'irga_warning_Tot(1)' %first cut if 1/2 have warning on then remove
    good = good & D(135,:) ~= 480; %'irga_warning_Tot(1)' Hard 480s - no sure but maybe calibration ?
    good = good & D(134,:) < 3600; %'csat_warning_Tot(1)' out of 7200
    %good = good & D(112,:) > minFlow(iSite);  %'irgaflow_Avg', l/min ==> need sufficient flow to reduce tube smearing
    good = good & D(110,:) > minPress(iSite);  %'irgapress_Avg(1)'
    good = good & D(110,:) < 120; %pressure should be in this range 
       
    %Use this later to remove bad FAST data - only sonic and irga based
    fast_good_irga_sonic=good;
    
    %NaN out covariances using sonic and irga
    D(63,~good)= NaN; %'Uz_co2_2(1)'
    D(64,~good)= NaN; %'Uz_h2o_2(1)'
    D(68,~good)= NaN; %'Ux_co2_2(1)'
    D(69,~good)= NaN; %'Ux_h2o_2(1)'
    D(72,~good)= NaN; %'Uy_co2_2(1)'
    D(73,~good)= NaN; %'Uy_h2o_2(1)'
    D(77,~good)= NaN; %'co2_Ts_2(1)' 
    D(79,~good)= NaN; %'h2o_Ts_2(1)'
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %bad irga
    %IRGA not working
    good = D(143,:) > -998; %'Red_alert_Tot(1)'
    good = good & D(144,:) > -998; %'Yellow_alert_Tot(1)'
    good = good & D(146,:) > 3600; % 'IRGA_on_Tot' irga is on for half the possible measurements
    good = good & D(135,:) < 3600; %'irga_warning_Tot(1)' %first cut if 1/2 have warning on then remove
    good = good & D(135,:) ~= 480; %'irga_warning_Tot(1)' Hard 480s - no sure but maybe calibration ?
    %good = good & D(112,:) > minFlow(iSite);  %'irgaflow_Avg', l/min ==> need sufficient flow to reduce tube smearing
    good = good & D(110,:) > minPress(iSite);  %'irgapress_Avg(1)'
    good = good & D(110,:) < 120; %pressure should be in this range 
    
    %Use this later to remove bad FAST data
    fast_good_irga = good;
    
    %NaN out irga
    D(75,~good)= NaN; %'co2_co2_2(1)'
    D(76,~good)= NaN; %'co2_h2o_2(1)'
    D(78,~good)= NaN; %'h2o_h2o_2(1)'
    D(85,~good)= NaN; %'co2_molar_Avg'
    D(86,~good)= NaN; %'h2o_molar_Avg'
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %SONIC is not working ==> Sonic fluxes are bad
    good = D(143,:) > -998; %'Red_alert_Tot(1)'
    good = good & D(134,:) < 3600; %'csat_warning_Tot(1)' out of 7200
    
    %should NaN out all SONIC based measurements
    D(60:74,~good)= NaN; %'Uz_Uz_2(1) through 'Uy_Ts_2(1)'
    D(77,~good)= NaN; %'co2_Ts_2(1)'
    D(79:84,~good)= NaN; %'h2o_Ts_2(1) through'Ts_1_Avg(1)'

    %Use this later to remove bad FAST data
    fast_good_sonic=good;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %GOES GOES GOES GOES GOES GOES GOES GOES GOES GOES GOES GOES GOES GOES GOES
    %(2)Alerts from goes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %'csat_warning and irga warnings - use this to remove -999s - important
    %for filtering
    
    bad = D(260,:) < -1; 
    D(260,bad)=NaN;
    
    bad = D(261,:) < -1; 
    D(261,bad)=NaN;
    
    %sonic and/or irga not working ==> sonic and irga measurements are bad
    good = D(269,:) > -998; %'GOES_red_alert'
    good = good & D(270,:) > -998; %'GOES_yellow_alert'
    good = good & D(272,:) > 3600; %GOES_irga_on
    good = good & D(261,:) < 3600; %'GOES_irga_warnings' %first cut if 1/2 have warning on then remove
    good = good & D(261,:) ~= 480; % calibration 
    good = good & D(260,:) < 3600; %'GOES_csat_warnings'
    %good = good & D(238,:) > minFlow(iSite); %'GOES_irgaflow_avg' l/min ==> need sufficient flow to reduce tube smearing
    good = good & D(236,:) > minPress(iSite); %'GOES_irgapress_avg'
    good = good & D(236,:) < 120; %'GOES_irgapress_avg'

    %NaN out covariances using sonic and irga
    D(189,~good)= NaN; %'GOES_Uz_co2'
    D(190,~good)= NaN; %'GOES_Uz_h2o'
    D(194,~good)= NaN; %'GOES_Ux_co2'
    D(195,~good)= NaN; %'GOES_Ux_h2o'
    D(198,~good)= NaN; %'GOES_Uy_co2'
    D(199,~good)= NaN; %'GOES_Uy_h2o'
    D(203,~good)= NaN; %'GOES_co2_Ts' 
    D(205,~good)= NaN; %'GOES_h2o_Ts'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %bad Irga
    good = D(269,:) > -998; %'GOES_red_alert'
    good = good & D(270,:) > -998; %'GOES_yellow_alert'
    good = good & D(272,:) > 3600; %GOES_irga_on
    good = good & D(261,:) < 3600; %'GOES_irga_warnings' %first cut if 1/2 have warning on then remove
    good = good & D(261,:) ~= 480; % calibration ?
    %good = good & D(238,:) > minFlow(iSite); %'GOES_irgaflow_avg' l/min ==> need sufficient flow to reduce tube smearing
    good = good & D(236,:) > minPress(iSite); %'GOES_irgapress_avg'
    good = good & D(236,:) < 120; %'GOES_irgapress_avg'

    %NaN out irga measurements
    D(201, ~good)= NaN; %'GOES_co2_co2'
    D(202, ~good)= NaN; %'GOES_co2_h2o'
    D(204, ~good)= NaN; %'GOES_h2o_h2o'
    D(211, ~good)= NaN; %'GOES_co2_avg'
    D(212, ~good)= NaN; %'GOES_h2o_avg'

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Now remove just bad sonic measurements
    good = D(269,:) > -998; %'GOES_red_alert'
    good = good & D(260,:) < 3600; %'GOES_csat_warnings'
    
    %NaN out all Sonic based measurements
    D(186:200, ~good)= NaN;
    D(203, ~good)= NaN;
    D(205:210, ~good)= NaN;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FAST %FAST  %FAST  %FAST  %FAST  %FAST  %FAST  %FAST  %FAST  %FAST 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %(3)FAST processed data
    %Assumption Alerts at data logger = alerts on fast processed data
    %==>NOTE: it is critical that the time stamp is being read correctly
    %fast data is already filtered in the fastflux program
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %bad irga and sonic
    D(30, ~fast_good_irga_sonic)= NaN;
    D(46:49, ~fast_good_irga_sonic)= NaN;
    
    %bad irga
    D(32:45, ~fast_good_irga)= NaN;

    %bad SONIC measurements
    D(5:29, ~fast_good_sonic)= NaN;
    
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
    n_sens_change = length(c);
    
    if n_sens_change == 2
        curr_multiplier=D(1,:)*NaN; 
        DateChangeR=datenum(DATE_INSTALLED(good,1))-towerYearStart(iSite);
        good = D(1,:) >= 0 & D(1,:) < DateChangeR(2,1);
        curr_multiplier(good)=c(1);
        good = D(1,:) >= DateChangeR(2,1) & D(1,:) <  max_D;
        curr_multiplier(good)=c(2);
    elseif n_sens_change == 1
        o=ones(1,size(D,2));
        curr_multiplier=o.*c(1);
    elseif n_sens_change > 2 || n_sens_change ==0
        disp('the code is not converting net radiation sensor voltages to physical units')
        pause
    end

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
    n_sens_change = length(c);
 
    if n_sens_change == 2
        curr_multiplier=D(1,:)*NaN; 
        good = D(1,:) >= 0 & D(1,:) < DateChangeR(2,1);
        curr_multiplier(good)=c(1);
        good = D(1,:) >= DateChangeR(2,1) & D(1,:) <  max_D;
        curr_multiplier(good)=c(2);
    elseif n_sens_change == 1
        o=ones(1,size(D,2)); 
        curr_multiplier=c(1).*o;
    elseif n_sens_change > 2 || n_sens_change ==0
        disp('the code is not converting net radiation sensor voltages to physical units')
        pause
    end

    %correct and replace data logger
    good = D(87,:) < 0;
    D(87,good)= D(87,good).* curr_multiplier(1,good);

    %correct and replace goes
    good = D(213,:) < 0;
    D(213,good)= D(213, good) .* curr_multiplier(1,good);

    
    %need this for adjustment
    Press=sitePress(iSite); %kPa - awf
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Section 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Site specific filtering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    if iSite == 1
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_DCBurnWorking(HEADER, D, Press, R_mol, Tc, 1);
    elseif iSite ==2 %DC Low Des
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_DC_LowDesWorking(HEADER, D, Press, R_mol, Tc, 1);
    elseif iSite ==3 % DC_Pinyon
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_PinyonWorking(HEADER, D, Press, R_mol, Tc, 1);
    elseif iSite ==4 %LR grassland
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_LRGrassWorking(HEADER, D, Press, R_mol, Tc, 1);
    elseif iSite ==5 %LR sage
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_LRSageWorking(HEADER, D, Press, R_mol, Tc, 1);
    elseif iSite ==6 %JamesRes
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_JamesWorking3(HEADER, D, Press, R_mol, Tc, 1);
    elseif iSite ==7
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_P301(HEADER, D, Press, R_mol, Tc, 1);
    elseif iSite ==8
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_SJER(HEADER, D, Press, R_mol, Tc, 1);
    elseif iSite ==9
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_Shorthair(HEADER, D, Press, R_mol, Tc, 1);
    elseif iSite ==10
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_Soaproot(HEADER, D, Press, R_mol, Tc, 1);
    else
        disp(['you have no site specific filters for ' siteName]) 
        disp('you may have introduced bad data into your output')
        pause
    end

    
    %----------------------------------------------------------------------
    %build up D so we can put in new calculations
    %----------------------------------------------------------------------
    column_D=length(D(1,:));
    row_D=length(D(:,1));
    Addrow_D=355-row_D;
    n=ones(Addrow_D, column_D).*NaN;
    D=[D;n]; 
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Section 4: Calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculation new variables 
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
    %horizontal wind vectors units = m/s
    %FAST
    u_mat=D(5,:);%'U mean, East to West'
    v_mat=D(6,:);%'V mean, North to South');
    
    %DATA LOGGER
    u_dl = D(81,:); %'Ux_1_Avg(1)');
    v_dl = D(82,:); %'Uy_1_Avg(1)');
    
    %goes
    u_goes = D(207,:); %'GOES_U_avg'; 
    v_goes = D(208,:); %'GOES_Uy_avg';
    
    %consensus
    u_consensus=uvw(1,:);
    v_consensus=uvw(2,:);
    
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
    
    %Consensus
    theta_consensus = (-atan2(v_consensus,u_consensus) + pi/2) *180/pi;
    theta_consensus(theta_consensus < 0)= theta_consensus(theta_consensus < 0)+360;
    
    theta_consensus = theta_consensus - 90 + SonicOrientation(iSite) ;
    theta_consensus(theta_consensus>360)=theta_consensus(theta_consensus>360)- 360;
    theta_consensus(theta_consensus < 0)= theta_consensus(theta_consensus < 0)+360;

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
    
    %GOES
    wspd_consensus = sqrt(u_consensus.^2 + v_consensus.^2);
    
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
    
    %(This should be T = Ts / (1 + 0.00032q ???? I changed this.  -awf)

    %% (1) could fill here first to use optimal data for calculation
    %% (2) H2o concentrations are from irga - which are not calibrated well
    
    %GET SOME DATA FOR THE CALCULATION
    %sonic T - Celcius
    Ts_mat = D(8,:);%'Ts');
    Ts_dl   = D(84,:);% 'Ts_1_Avg(1)';
    Ts_goes = D(210 ,:);%'GOES_Ts_avg');
    
    %make sure Ts is in deg C
    if mednn(Ts_mat) > 100
        Ts_mat = Ts_mat - Tc;
    end

    if mednn(Ts_dl) > 100
        Ts_mat=Ts_dl-Tc;
    end

    if mednn(Ts_goes)>100
        Ts_goes=Ts_goes-Tc;
    end
    
    %!!!!!!!!!!!!!!!!CHECK THE UNITS!!!!!!!!!!!!!!!!!
    %IRGA water vapor concentrations
    %These are water mixing ratios mmol of h2o/mol dry air:
    h2o_mat  = D(42 ,:); %'H2O Closed Path MEAN (mmol/mol)');
    h2o_dl   = D(86 ,:); %'h2o_molar_Avg');
    h2o_goes = D(212 ,:); %'GOES_h2o_avg');
    %h2o_consensus=h2o_consensus; this is mostly the hmp filled with irga
    
    %If you want mole fraction mmol of h2o/mol of wet air
    molefrach2o_mat  = D(42 ,:)./(1+(D(42 ,:)/1000)); 
    molefrach2o_dl   = D(86 ,:)./(1+(D(86 ,:)/1000));
    molefrach2o_goes = D(212 ,:)./(1+(D(212 ,:)/1000));
    molefrach2o_consensus = h2o_consensus./(1+(h2o_consensus(1,:)/1000));
    
    %Kaimal and Gaynor, 1991 - should be mole fraction not mixing ratio - awf
    %FAST
    Tkact_mat   = (Ts_mat + Tc)  ./ (1 + 0.00032 * molefrach2o_consensus); %[K]
    
    %Add to D and Header
    D(330,:)=Tkact_mat; %K
    HMERGE(330,1)={'Tkact_mat'};
    
    
    %DATALOGGER
    Tkact_dl    = (Ts_dl + Tc)   ./ (1 + 0.00032 * molefrach2o_consensus);%K
    
    %Add to D and Header
    D(331,:)=Tkact_dl;
    HMERGE(331,1)={'Tkact_dl'};
    
    %GOES
    Tkact_goes  = (Ts_goes + Tc) ./ (1 + 0.00032 * molefrach2o_consensus);
    
    %Add to D and Header
    D(332,:)=Tkact_goes;
    HMERGE(332,1)={'Tkact_goes'};
    
    %----------------------------------------------------------------
    %% Moist Air Density using dried Sonic Temperature
    %----------------------------------------------------------------
    %This is the best air density to use because it is at the sonic - where
    %the sampling is taking place. Temp unnits are in Kelvin - awf
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
    

    %use this for calculations
    rhom_consensus = Press ./ (R_mol * (T_consensus+Tc))  * 1e3; %units mol moist air per m3 
    
    %----------------------------------------------------------------
    %% Partial pressure of water vapor - same as e=Press*(h2o/1000) -awf
    %----------------------------------------------------------------
    %should use mole fraction for h2o - awf
    Pw_mat  = rhom_ma_mat  .* molefrach2o_mat  .* Tkact_mat  * R_mol * 1e-6; %kPa
    Pw_dl   = rhom_ma_dl   .* molefrach2o_dl   .* Tkact_dl   * R_mol * 1e-6; %kPa
    Pw_goes = rhom_ma_goes .* molefrach2o_goes .* Tkact_goes * R_mol * 1e-6; %kPa

    Pw_consensus = rhom_consensus .* molefrach2o_consensus .* (T_consensus+Tc) * R_mol * 1e-6; %kPa
    
    %----------------------------------------------------------------
    %% Partial pressure of Dry Air - Dalton's law of partial pressure
    %---------------------------------------------------------------- 
    Pa_mat  = Press - Pw_mat; %kPa
    Pa_dl   = Press - Pw_dl;
    Pa_goes = Press - Pw_goes;

    Pa_consensus = Press - Pw_consensus;
    
    %----------------------------------------------------------------
    %% Density of Dry Air
    %----------------------------------------------------------------
    rhom_da_mat  = 1e3 * Pa_mat  ./ (R_mol *  Tkact_mat); %mol/m^3
    rhom_da_dl   = 1e3 * Pa_dl   ./ (R_mol *  Tkact_dl);
    rhom_da_goes = 1e3 * Pa_goes ./ (R_mol *  Tkact_goes);

    rhom_da_consensus = 1e3 * Pa_consensus ./ (R_mol *  (T_consensus+Tc));
    
    %The fast calculated fluxes use bad temp and water vapor to calculate
    %dry air density, which is used to calculate fluxes
    
    %"need_rho_fix" is a logical row vector created in the site specfic
    %calculations that tells us if the fast data air temperature or water
    %vapor was bad (air temperature or water vapor were used to calculate dry air density)
    %==>we can adjust these fluxes
    
    %This is the best consensus dry air density/dry air density used in the
    %fast calculation
    %use this to adjust matlab processed fluxes
    densityfix=rhom_da_consensus./D(48,:);
    
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
    %% (Water vapor molar density? - awf) Wet air molar density (mol/m3) and total mass density (kg/m3)
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
    %P varies with altitude - awf - 2012
    
    % get sensible heat flux from goes data
    % use the equation:    H = m_a * rho_a * CpAirc * <w'.T'>
    % where m_a = molar mass of air 28.966/1000
    % rho_a = density of air (mol m-3) calcd by rho = P/RT assum P is const (101.3 kPa) and R = 8.314 J/K/mol
    % CpAirc is the heat content of dry air (1004.67 J kg-1 K-1)
    % <w'.T'> is the 30 minute average of the vertical wind and the sonic
    % temperature
    
    %dry air density adjustment to mat proc
    D(29,need_rho_fix) = densityfix(1,need_rho_fix).*D(29,need_rho_fix);
    
    %----------------------------------------------------------------
    % declare some constants initially
    m_a= Mw_da / 1000; %units of kg dry air / mol
    %----------------------------------------------------------------
    
    %DATA LOGGER
    sh_dl   = m_a * CpAirc * rhom_da_consensus.* D(65,:); %units of W m-2; %D(65.:)='Uz_Ts_2(1)';

    %Add to D and Header
    D(339,:)= sh_dl;
    HMERGE(339,1)={'sh_dl W m-2'};
    
    
    %GOES
    sh_goes = m_a * CpAirc * rhom_da_consensus .* D(191,:); %units of W m-2;

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
   
    %dry air density adjustment to mat proc
    D(30,need_rho_fix) = densityfix(1,need_rho_fix).*D(30,need_rho_fix);
    
    %Ts should be in Celcius
    %LOGGER
    %calculate Latent heat of vaporization using the sonic temperature 
    Lv_dl   = (2.501-0.00237.*(T_consensus))*10^6; %units: J/kg
    
    %calculate latent heat
    lh_dl   = 1e-6 * 18 * rhom_da_consensus.*Lv_dl .* D(64,:); %'Uz_h2o_2(1)'
    
    %Add to D and Header
    D(343,:)=lh_dl;
    HMERGE(343,1)={'lh_dl'};
    
    %GOES
    %calculate Latent heat of vaporization using the sonic temperature
    Lv_goes = (2.501-0.00237.*(T_consensus))*10^6; %units: J/kg
    
    %calculate latent heat
    lh_goes = 1e-6 * 18 * rhom_da_consensus.*Lv_goes .* D(190,:); %'GOES_Uz_h2o'

    %Add to D and Header
    D(344,:)=lh_goes;
    HMERGE(344,1)={'lh_goes'};

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
    
    %dry air density adjustment to mat proc
    D(46,need_rho_fix) = densityfix(1,need_rho_fix).*D(46,need_rho_fix);
    
    %datalogger
    Fco2_dl = rhom_da_consensus.* D(63,:);%'Uz_co2_2(1)'; 
    
    %save the data
    D(345,:)=Fco2_dl;
    HMERGE(345,1)={'Fco2_dl'};
    
    %goes
    Fco2_goes  = rhom_da_consensus.* D(189,:); %'GOES_Uz_co2';
    
    %save the data
    D(346,:)=Fco2_goes;
    HMERGE(346,1)={'Fco2_goes'};
    
    %----------------------------------------------------------------
    %% H2O fluxes
    %----------------------------------------------------------------
    
    %dry air density adjustment to mat proc
    D(47,need_rho_fix) = densityfix(1,need_rho_fix).*D(47,need_rho_fix);
    
    %data logger
    Fh2o_dl = rhom_da_consensus.* D(64,:);%'Uz_h2o_2(1)' 
    
    %save the data
    D(347,:)=Fh2o_dl;
    HMERGE(347,1)={'Fh2o_dl'};
    
    %goes
    Fh2o_goes = rhom_da_consensus.* D(190,:);%'GOES_Uz_h2o'
    
    %save the data
    D(348,:)=Fh2o_goes;
    HMERGE(348,1)={'Fh2o_goes'};
    %now correct for the delay and smearing in the datalogger and goes data using the Shaw Tau model

    %     fh2o_dl   = fh2o_dl.*gain_fh2o;
    %     fh2o_goes = fh2o_goes.*gain_fh2o;
    
    %----------------------------------------------------------------
    %Add consensus calculations to D
    %----------------------------------------------------------------
    D(349,:)=T_consensus;
    HMERGE(349,1)={'T_consensus'};
    
    D(350,:)=h2o_consensus;
    HMERGE(350,1)={'h2o_consensus'};
    
    D(351,:)=need_rho_fix;
    HMERGE(351,1)={'need_rho_fix'};
    
    D(352,:)=rhom_consensus;
    HMERGE(352,1)={'rhom_consensus'};
    
    D(353,:)=rhom_da_consensus;
    HMERGE(353,1)={'rhom_da_consensus'};
    
    D(354,:)=theta_consensus;
    HMERGE(354,1)={'theta_consensus'};
    
    D(355,:)=wspd_consensus;
    HMERGE(355,1)={'wspd_consensus'};
    
    %----------------------------------------------------------------------
    %redo HEADER because we added new rows to the data D
    %----------------------------------------------------------------------
    HEADER = char(HMERGE);
    %----------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Section 5 : Mop-up the outliers that occurred in calculations 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    if iSite == 1
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_DCBurnWorking(HEADER, D, Press, R_mol, Tc, 2);
    elseif iSite ==2
        [D, T_consensus, h2o_consensus, need_rho_fix] = Site_specific_DC_LowDesWorking(HEADER, D, Press, R_mol, Tc, 2);
    elseif iSite ==3
        [D, T_consensus, h2o_consensus, need_rho_fix] = Site_specific_PinyonWorking(HEADER, D, Press, R_mol, Tc, 2);
    elseif iSite ==4 %LR grassland
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_LRGrassWorking(HEADER, D, Press, R_mol, Tc, 2);
    elseif iSite ==5
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_LRSageWorking(HEADER, D, Press, R_mol, Tc, 2);
    elseif iSite ==6
        [D, T_consensus, h2o_consensus, need_rho_fix] =Site_specific_JamesWorking3(HEADER, D, Press, R_mol, Tc, 2);
    elseif iSite ==7
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_P301(HEADER, D, Press, R_mol, Tc, 2);
    elseif iSite ==8
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_SJER(HEADER, D, Press, R_mol, Tc, 2);
    elseif iSite ==9
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_Shorthair(HEADER, D, Press, R_mol, Tc, 2);
    elseif iSite ==10
        [D, T_consensus, h2o_consensus, need_rho_fix, uvw] = Site_specific_Soaproot(HEADER, D, Press, R_mol, Tc, 2);
    else
        disp(['you have no site specific filters for ' siteName]) 
        disp('you may have introduced bad data into your output')
        pause
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Section 6    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Here is where we combined
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %Each dataset has been filtered and corrected for sensor calibrations 
    %Now we will combined these different data streams into 1 data table
    %We combined data streams in the following order: 1) fast, 2)datalogger, 3)goes, 4) consensus - merged above, 5) soil logger 
    
    %What is the variable? orginized by column
    Cleaned_Header={};
    
    %Here is the data for the variable in Cleaned_Header
    ncols=size(D,2);
    Cleaned_D=ones(29, ncols)*NaN;
    
    %What data stream did the data come from?
    DataSource=ones(29, ncols)*NaN; %1 = processed; 2=datalogger; 3=GOES
    
    %Do things look OK?
    if combinInfo ==1
        [N_obs, headerN_obs, D_obs_out] = combinedinfo(siteName, D, D_obs_out, N_obs, 2);
    end
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
    %homogenize
    disp('homogenize Ustar with robustfit')
    isn=isnan(D(21,:)) ==1 | isnan(D(322,:)) == 1;
    %P=polyfit(D(21,~isn),D(322,~isn),1);
    %Dadj=(D(322,:)-P(1,2))./P(1,1);
    B = robustfit(D(21,~isn),D(322,~isn)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(322,:)-B(1,1))./B(2,1);
    D(322,:)=Dadj;
    
    isn=isnan(D(21,:)) ==1 | isnan(D(323,:)) == 1;
    %P=polyfit(D(21,~isn),D(323,~isn),1);
    %Dadj=(D(323,:)-P(1,2))./P(1,1);
    B = robustfit(D(21,~isn),D(323,~isn)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(323,:)-B(1,1))./B(2,1);
    D(323,:)=Dadj;
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
    %homogenize
    isn=isnan(D(327,:)) ==1 | isnan(D(328,:)) == 1;
    P=polyfit(D(327,~isn),D(328,~isn),1);
    Dadj=(D(328,:)-P(1,2))./P(1,1);
    D(328,:)=Dadj;
    
    isn=isnan(D(327,:)) ==1 | isnan(D(329,:)) == 1;
    P=polyfit(D(327,~isn),D(329,~isn),1);
    Dadj=(D(329,:)-P(1,2))./P(1,1);
    D(329,:)=Dadj;
    
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
    %homogenize
    isn=isnan(D(336,:)) ==1 | isnan(D(337,:)) == 1;
    P=polyfit(D(336,~isn),D(337,~isn),1);
    Dadj=(D(337,:)-P(1,2))./P(1,1);
    D(337,:)=Dadj;
    
    isn=isnan(D(336,:)) ==1 | isnan(D(338,:)) == 1;
    P=polyfit(D(336,~isn),D(338,~isn),1);
    Dadj=(D(338,:)-P(1,2))./P(1,1);
    D(338,:)=Dadj;
    %****************************************
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
    %homogenize
    isn=isnan(D(333,:)) ==1 | isnan(D(334,:)) == 1;
    P=polyfit(D(333,~isn),D(334,~isn),1);
    Dadj=(D(334,:)-P(1,2))./P(1,1);
    D(334,:)=Dadj;
    
    isn=isnan(D(333,:)) ==1 | isnan(D(335,:)) == 1;
    P=polyfit(D(333,~isn),D(335,~isn),1);
    Dadj=(D(335,:)-P(1,2))./P(1,1);
    D(335,:)=Dadj;
    %****************************************
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
    %homogenize
    isn=isnan(D(29,:)) ==1 | isnan(D(339,:)) == 1;
    P=polyfit(D(29,~isn),D(339,~isn),1);
    Dadj=(D(339,:)-P(1,2))./P(1,1);
    D(339,:)=Dadj;
    
    isn=isnan(D(29,:)) ==1 | isnan(D(340,:)) == 1;
    P=polyfit(D(29,~isn),D(340,~isn),1);
    Dadj=(D(340,:)-P(1,2))./P(1,1);
    D(340,:)=Dadj;
    %****************************************
    
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
    %homogenize
    isn=isnan(D(30,:)) ==1 | isnan(D(343,:)) == 1;
    %P=polyfit(D(30,~isn),D(343,~isn),1);
    %Dadj=(D(343,:)-P(1,2))./P(1,1);
    disp('using robustfit to homogenize LE')
    B = robustfit(D(30,~isn),D(343,~isn)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(343,:)-B(1,1))./B(2,1);
    D(343,:)=Dadj;
    
    isn=isnan(D(30,:)) ==1 | isnan(D(344,:)) == 1;
    %P=polyfit(D(30,~isn),D(344,~isn),1);
    %Dadj=(D(344,:)-P(1,2))./P(1,1);
    B = robustfit(D(30,~isn),D(344,~isn)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(344,:)-B(1,1))./B(2,1);
    D(344,:)=Dadj;
    %****************************************
    
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
    %homogenize
    isn=isnan(D(35,:)) ==1 | isnan(D(85,:)) == 1;
    P=polyfit(D(35,~isn),D(85,~isn),1);
    Dadj=(D(85,:)-P(1,2))./P(1,1);
    D(85,:)=Dadj;
    
    isn=isnan(D(35,:)) ==1 | isnan(D(211,:)) == 1;
    P=polyfit(D(35,~isn),D(211,~isn),1);
    Dadj=(D(211,:)-P(1,2))./P(1,1);
    D(211,:)=Dadj;
    %****************************************
    
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
    %12. H2O concentration - mixing ratio
    %----------------------------------------------------------------------
    %homogenize
    isn=isnan(D(42,:)) ==1 | isnan(D(86,:)) == 1;
    P=polyfit(D(42,~isn),D(86,~isn),1);
    Dadj=(D(86,:)-P(1,2))./P(1,1);
    D(86,:)=Dadj;
    
    isn=isnan(D(42,:)) ==1 | isnan(D(212,:)) == 1;
    P=polyfit(D(42,~isn),D(212,~isn),1);
    Dadj=(D(212,:)-P(1,2))./P(1,1);
    D(212,:)=Dadj;
    %****************************************
    
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
    %split into neg fluxes and pos fluxes before adjusting the gain
    
    %we will use just ust filtered observations 
    turb = Cleaned_D(4,:)>0.3;
    
    %homogenize
    Dadj(1,:)=D(345,:).*NaN;
    isn=isnan(D(46,:)) ==1 | isnan(D(345,:)) == 1;
    
    reg = D(345,:) <= 0 & D(46,:) <= 0 & turb ==1 & ~isn;
    %P=polyfit(D(46,reg),D(345,reg),1);
    %least squares regression with intercept through zero
    A=D(46,reg)';
    B=D(345,reg)';
    P=A\B; % matlab gives you the least squares solution to the A P = B
    neg = D(345,:) <= 0;
    Dadj(1,neg)=D(345,neg)/P(1,1);
    
    reg = D(345,:) > 0 & D(46,:) > 0 & turb ==1 & ~isn;
    %P=polyfit(D(46,reg),D(345,reg),1);
    A=D(46,reg)';
    B=D(345,reg)';
    P=A\B; % matlab gives you the least squares solution to the A P = B
    pos = D(345,:) > 0;
    Dadj(1,pos)=D(345,pos)./P(1,1);
    
    D(345,:)=Dadj;
    
    Dadj(1,:)=D(346,:).*NaN;
    isn=isnan(D(46,:)) ==1 | isnan(D(346,:)) == 1;
    reg = D(346,:) <= 0 & D(46,:) <= 0 & turb ==1  & ~isn;
    %P=polyfit(D(46,reg),D(346,reg),1);
    A=D(46,reg)';
    B=D(346,reg)';
    P=A\B; % matlab gives you the least squares solution to the A P = B
    neg = D(346,:) <= 0;
    Dadj(1,neg)=D(346,neg)/P(1,1);
    
    reg = D(346,:) > 0 & D(46,:) > 0 & turb ==1  & ~isn;
    %P=polyfit(D(46,reg),D(346,reg),1);
    A=D(46,reg)';
    B=D(346,reg)';
    P=A\B; % matlab gives you the least squares solution to the A P = B
    pos = D(346,:) > 0;
    Dadj(1,pos)=D(346,pos)/P(1,1);
    D(346,:)=Dadj;
    %****************************************
    
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
    %we will use just ust filtered observations 
    turb = Cleaned_D(4,:)>0.3;
    
    %homogenize
    isn=isnan(D(47,:)) ==1 & turb ==1 | isnan(D(347,:)) == 1 & turb ==1;
    %P=polyfit(D(47,~isn),D(347,~isn),1);
    %Dadj=(D(347,:)-P(1,2))./P(1,1);
    disp('using robustfit to homogenize Fh2o')
    B = robustfit(D(47,~isn),D(347,~isn)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(347,:)-B(1,1))./B(2,1);
    D(347,:)=Dadj;
    
    isn=isnan(D(47,:)) ==1 & turb ==1 | isnan(D(348,:)) == 1 & turb ==1;
    %P=polyfit(D(47,~isn),D(348,~isn),1);
    %Dadj=(D(348,:)-P(1,2))./P(1,1);
    disp('using robustfit to homogenize Fh2o')
    B = robustfit(D(47,~isn),D(348,~isn)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(348,:)-B(1,1))./B(2,1);
    D(348,:)=Dadj;
    %****************************************
    
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
    Cleaned_D(35:44,:)=D(282:291,:);%dl
    Cleaned_D(45:53,:)=D(293:301,:);%dl
    Cleaned_D(54:59,:)=D(314:319,:);%dl
    Cleaned_D(60:61,:)=D(320:321,:);%dl
    
    %Where is the data from?
    DataSource(35:61,:)=5;

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
    Cleaned_Header(54,1)={'DelT_C(1)'}; 
    Cleaned_Header(55,1)={'DelT_C(2)'}; 
    Cleaned_Header(56,1)={'DelT_C(3)'}; 
    Cleaned_Header(57,1)={'DelT_C(4)'}; 
    Cleaned_Header(58,1)={'DelT_C(5)'}; 
    Cleaned_Header(59,1)={'DelT_C(6)'}; 
    Cleaned_Header(60,1)={'AirT'}; 
    Cleaned_Header(61,1)={'Snow Depth'}; 
    
    %Flip to maintain format
    Cleaned_D=Cleaned_D';
    DataSource=DataSource';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Section 7
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Save
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    fnamecleaned = [combineRootDir  siteName '_cleaned'];

    if saveon ==1
        disp('You saved your results: check')
        fnamecleaned
        save(fnamecleaned, 'Cleaned_D', 'Cleaned_Header', 'DataSource', 'D', 'HMERGE', 'Day');
    else
        disp('you did not save your results')
    end

end