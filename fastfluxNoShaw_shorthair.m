function fastfluxNoShaw_shorthair(Sites2Proc)

    
% current version 2/3/11 aek
% simplified, reordered, and changed Shaw Tau to covarmax 5/14/12 awf
%% Tasks
% - add error messages to RunInfo, esp. from subroutines
% - fix psd, csd functions in taucalc.m
% - get macros to read "data2" columns
% - fix column headings for "_avg" and " avg" column header names

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           ~~~~  Part 1. Initialilizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except Sites2Proc;
path(path, 'C:\towerData\ProcessingScripts\subroutines');
global sites sitePress towerStartDates iSite towerYearStart VIsites MVL_Universal Tc R_mol
% global siteAlt
global hostname dirSep 
global outputFileRootName
global fastRootDir inputRootDir
global inputRootDirList
%% flux processing parameters
global procInt qTaus tauRefc tauRefh Delay_estimate Delay_estimate_co2 Delay_estimate_h2o Tau_estimate_co2 Tau_estimate_h2o
%%
var_defs();
Day = date;
%%
close all
tic

%if nargin==0
%    Sites2Proc = 1:length(sites);
%end

diary_filename = [fastRootDir, 'fast_log_', Day, char(sites(Sites2Proc(1)))];
diary(diary_filename);
format long g


%           ~~~~  Start the Site Loop
for iSite = Sites2Proc      %L1
 
    Site = char(sites(iSite));
    JDStart = towerStartDates(iSite);
    disp(['Todays date is ' Day]);
    disp(['Getting Ready to Process ' Site]);

    close all
    clear D HEADER;
    clear UVWTMEAN UVWTVAR UVWTSKEW UVWTKURT THETA COVUVWT USTAR HBUOYANT;
    clear H30MIN D30_NEW D30MIN JDInterval TIME HHMM TAUC TAUH DELAYC DELAYH;
    clear d30 CO2 FCO2 H2O FH2O RHOM HLATENT HSENSIBLE GAINC GAINH;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           ~~~~  End Part 1. Initialilizations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Do you want to search for a new Tau and Delay - awf
    TauandDelay = 0; %1 = yes; 0 = no
    %using this to store Tau estimates
    TCR=[];
    
    
    %These are the Taus and delays you will use.  They are site specific.  These are set in var_def
    %- awf 5/31/2012
    delayRefc = Delay_estimate_co2(iSite);    
    delayRefh = Delay_estimate_h2o(iSite);

    tauRefc = Tau_estimate_co2(iSite);
    tauRefh = Tau_estimate_h2o(iSite);
    
    %removed Shaw 5/31/2012 - awf
       
    %*******************************************************
    %               Part 2.  BEGINNING OF USER INPUTS
    %*******************************************************

    %%% fix the processing constants
    RunInfo = [date,...
            ['California Data ' Site],...
            ['Processed on ' hostname],...
            'Estimating delay based on max covariance between w and scalar',...
            'Using mean delay of Delay(co2) = delayRefc and delay(h2o) = delayRefh (see var_def file)',...
            'Eventually, we will tightened window around these delays for looking for maximum covariance',...
            'to 2 seconds about these delays.'];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINE INPUT/OUTPUT FILES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fastSiteDir = [fastRootDir Site dirSep];
    OutputFileDir    = [fastSiteDir 'proc' dirSep];
    if exist(OutputFileDir,'dir') ~= 7
        mkdir (OutputFileDir);
    end
    
    %*******************************************************
    %               End Part 2.  END OF USER INPUTS
    %*******************************************************   
    %*******************************************************
    %     Part 3.  SET UP FILE STRUCTURE AND OUTPUT FILES
    %*******************************************************
    % output files
    fileout    = [OutputFileDir outputFileRootName Site '.mat'];
    fileoutbak = [OutputFileDir outputFileRootName Site '_' Day '.mat'];
    fileout_one_array_old    = [OutputFileDir outputFileRootName Site '_one_array_old_' Day '.mat'];  %save old version of one array
    fileout_one_array    = [OutputFileDir outputFileRootName Site '_one_array.mat'];   %new version of one array (joined old and new)
    fileout_one_arraybak = [OutputFileDir outputFileRootName Site '_' Day '_one_array.mat'];  %back up version of one array (joined old and new)

    % Write output to an ASCII file that can be read by Excel
    fileoutascii   = [OutputFileDir outputFileRootName Site '.asc'];
    channelfile    = [OutputFileDir outputFileRootName Site '.channels'];

    %%% add to RunInfo
    RunInfo=[RunInfo,['Input Directory: ' fastSiteDir],['Output .mat file: ' fileout],['Output .mat file backup: ' fileoutbak],['Output file, one array: ' fileout_one_array],['Output file, one array backup: ' fileout_one_arraybak],['Output ASCII file: ' fileoutascii],['List of Channels: ' channelfile]]; %#ok<*AGROW>

    %%% This is used for plotting the results of taucalc.m
    PlotInfo.figure  = 1;
    PlotInfo.sub1    = '(2,1,1)';
    PlotInfo.sub2    = '(2,1,2)';
    PlotInfo.title1   = 'Co2 Closed Path:  ';
    PlotInfo.title2   = 'Co2: Ts as Reference';
    PlotInfo.legend  = {'Measured Attenuated Signal','Reference Signal','Degraded Reference signal'};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONTINUE JOB OR START NEW JOB
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% First Deal with the the One Array data
    % 1. Copy last version of FluxesSITE_one_array to
    % FluxesSITE_one_array_old_date

    if exist(fileout_one_array,'file') == 2
        copyfile(fileout_one_array,fileout_one_array_old);
        eval(['load ' fileout_one_array_old ';']);

        %2. Now read this file
        DOLD = D;
        %HOLD = HEADER;
        HOLD =MVL_Universal;
        orig_length = size(D,2);
        clear D HEADER;

        %2b Check for the possiblity that the time fields have been updated but the
        % data fields are all full of -999.
        first_data_column = DOLD(4,:); 
        last_good_col = find(first_data_column~=-999, 1, 'last' );
        block = DOLD(4:end,last_good_col+1:end);
        if isempty(find(block~=-999, 1))
            DOLD = DOLD(:,1:last_good_col);
        end

        ncols_removed = orig_length-size(DOLD,2); %awf
        disp(['found and removed ' num2str(ncols_removed) ' trailing cols of valid times but -999 fills in data fields']);

        %3 Identify 2nd row of data as PREVTIMES and find out last non-Nan value      
        PREVTIMES = DOLD(2,:);
        last_good_ix = find(PREVTIMES>-999 & ~isnan(PREVTIMES), 1, 'last' ); 
        last_good_time = PREVTIMES(last_good_ix);
    else         
        last_good_time = [];
    end 
    
    %4. Now Go to the raw data files and work out when the new interval begins
    %and ends
    ts_data_Dir = [fastSiteDir 'ts_data' dirSep];

    ts_data_Dir_List = dir(ts_data_Dir);

    %start with ts_data....go through names
    numbergoodfiles=0;
    for ifile = 1:length(ts_data_Dir_List)
        cname = ts_data_Dir_List(ifile).name;
        if length(cname)==27 && strcmpi(cname(17:27),'ts_data.tob')
            numbergoodfiles = numbergoodfiles+1;
            YYYY = str2num(cname(1:4)); %#ok<ST2NM>
            MM = str2num(cname(6:7)); %#ok<ST2NM>
            DD = str2num(cname(9:10)); %#ok<ST2NM>
            daysAvail(numbergoodfiles) = datenum(YYYY,MM,DD) - towerYearStart(iSite);
        end
    end
    
    daysAvail=sort(daysAvail);
    firstDayAvail = daysAvail(1);
    lastDayAvail = daysAvail(end);

    if ~isempty(last_good_time)   % if no data processed yet
        disp(last_good_time)
        firstDayNotProc = daysAvail( find( daysAvail > (last_good_time-1), 1) ); %subtract 1 so that a partially processed file will now be fully processed
    else
        firstDayNotProc = firstDayAvail;
    end
    
    %changed this on 7/5/2012 to get the whole last day -awf
    numfilestoprocess = lastDayAvail - firstDayNotProc + 1;
    
    %5. Create i2Proc
    JDEnd = lastDayAvail+1;
    JDInterval = JDStart : procInt : JDEnd;
    i2Proc = find(JDInterval>=firstDayNotProc & JDInterval<=lastDayAvail);  % units: index to JDInterval
    int2Proc = JDInterval(i2Proc);   %New Variable Added Units Time in ED
    
    %6 Print Some Information to screen
    disp(['Total Number of files in ts_data directory = ' num2str(numbergoodfiles)]);
    disp(['Last Day processed was ' num2str(firstDayNotProc)  ' (' cal_ed2date(firstDayNotProc)   ') '  ]);
    disp(['About to process days  ' num2str(firstDayNotProc) ' ('  cal_ed2date(firstDayNotProc) ') to ' num2str(lastDayAvail)  ' ('  cal_ed2date(lastDayAvail) ')' ]);
    disp(['Number of files to process ' num2str(numfilestoprocess)]);
    disp(['Number of 30 Minute Intervals to process ' num2str(length(i2Proc))]);

    
    %--------------------------------------------------------------------------
    %% 7 Derive the variable names from D and add arrays of -999 to extend them
    %--------------------------------------------------------------------------
    % 7 (A) Deal with the fast data first. There are ~58 Rows of fast data
    % (including the Counter) and about ~23 variable groups. In the one_array
    % matrix DOLD, we must assign the variable groups back to the individual
    % array rows. This is done below. We deal with the slow (30 minute data
    % separately in 7(B) )

    %We already read in old header (HOLD)- it should match the universal
    %header once the code has run though ==> may delete this section in
    %future runs - awf
    
    %We want a dataset that has the same rows and headers
    
    ALLVARS = upper(MVL_Universal);%MVL_Universal is from var_def file
    nrows2add = length(int2Proc);
    
    %SAVE JUST SLOW 
    SLOW_D30= [];
    
    % dataset to the universal header
    D=[];
    nrows=length(ALLVARS);
        
    if exist(fileout_one_array,'file') == 2
        ncol=nrows2add;
        Dext=ones(nrows,ncol).*NaN;
        D=[DOLD, Dext];      
    else
        ncol=nrows2add;
        D=ones(nrows,ncol).*NaN; %if D does not exist ==> we need to make D - it should be nrows2add wide and length(ALLVARS) high
    end


    %used to figure out how much data we are going to redo
    DOLD_exist=exist('DOLD', 'var');
    
    if DOLD_exist == 0;
        nrows_DOLD = 0;
        orig_length = 0;
        ncols_removed =0;
    else
        nrows_DOLD = size(DOLD,1);
    end
    
    %size of new data set
    Dsize=size(D);
    
    disp([ ' Size of Old Data array was ' num2str(nrows_DOLD) ' x ' num2str(orig_length)]);
    disp([ ' Size of Old Data array after blank trailing columns removed was ' num2str(nrows_DOLD) ' x ' num2str(orig_length-ncols_removed)]);
    disp([ ' Size of New Data array is ' num2str(Dsize(1)) ' x ' num2str(Dsize(2))]);
      
        
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('%%%%    Now begin Processing intervals                       %%%%%%%');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% BEGIN PROCESSING THE INTERVALS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for ii=i2Proc
        D(1,ii)=ii; % a counter
        current_time = JDInterval(ii);
        TIME = current_time + 1/48; % points to the end of the 30 min interval
        D(2,ii)= TIME;
        hh = floor((rem(TIME, 1)*24)+(10/86400));
        mm = floor((rem(TIME, 1)*24*60)+(10/86400)) - hh*60; % we had rounding errors so I added 10 seconds to get floor to work - awf
        HHMM =  hh*100+mm;
        D(3,ii)= HHMM;
        disp(num2str([ ii HHMM JDInterval(ii) JDInterval(ii)+ procInt]));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Now, Read the data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %d=[];h=[];dflux=[];hflux=[];dvariance=[];hvariance=[];
        %it=[];ic=[];ih=[];SONDIAG=[];IRGADIAG=[];SF=[]; Commented out on
        %08/18/09 by RGA


        % Read in 40 minutes - 5 before the half hour and 5 after -
        
        % looks like 1 min to me???? - awf??
        
        % these will be truncated to calculate 30 minute fluxes, but
        % need the to get rid of filter edge effects for the closed
        % path.
        % Need to get the file list here for the interval to process

        disp([JDInterval(ii)-1/1440  JDInterval(ii)+procInt+1/1440 300]);
        [h,d]                 = concat5000([fastSiteDir 'ts_data' dirSep  ],[],[JDInterval(ii)-1/1440  JDInterval(ii)+procInt+1/1440 300]);
        [hflux,dflux]         = concat5000([fastSiteDir 'flux' dirSep     ],[],[JDInterval(ii)+procInt    JDInterval(ii)+procInt+1/24/3600 900]);
        [hvariance,dvariance] = concat5000([fastSiteDir 'variance' dirSep ],[],[JDInterval(ii)+procInt    JDInterval(ii)+procInt+1/24/3600 900]);
       
        if isempty(d)
            disp('ts data not being read.  Check concat5000 line in FSC program.')
            %pause;
        end
        if isempty(dvariance)
            disp('variance data not being read.  Check concat5000 line in FSC program.')
            %pause;
        end
        if isempty(hflux)
            disp('flux data not being read.  Check concat5000 line in FSC program.')
            %pause;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PROCESS THE FAST DATA HERE   hg(1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Begin with a requirement for 20 minutes of data

        if ~isempty(d)
            d(isnan(d))=-999;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ROTATE WINDS/CALCULATE MOMENTUM FLUX
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            s  = chset(d,h,'SECONDS');
            ns = chset(d,h,'NANOSECONDS');

            DayShift = towerYearStart(iSite) - datenum(1990,1,0);
    
            time_sec = s+ns/1e9; % turn to seconds
            fast_sec=max(time_sec)-60; %Line 433 accoutns for a 1 min edge on the interval - subtract this 1 min edge to get to time for the interval
           
            %This should agree with the data logger assuming that the data logger places time stamp on the end of the
            %30 min period in units of seconds from datenum(1990,1,0). -awf
            time_day = time_sec/24/3600-DayShift;  

            % calculate sampling frequency
            dt = median( diff(time_sec) ) ;
            SF = 1/dt;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% parameters for calculating time constants and delays
            PFlux.sf = SF;
            ParamsCo2.sf = SF;  ParamsH2o.sf = SF;
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            uvwt = [chset(d,h,'Ux_1(1)');chset(d,h,'Uy_1(1)');chset(d,h,'Uz_1(1)');chset(d,h,'Ts_1(1)')];
            Ts   =  chset(d,h,'Ts_1(1)')+Tc;
            diag =  chset(d,h,'diag_1(1)');

            %%% CR 5000 version for P301
%             uvwt = [chset(d,h,'Ux_1');chset(d,h,'Uy_1');chset(d,h,'Uz_1');chset(d,h,'Ts_1')];
%             Ts   =  chset(d,h,'Ts_1')+Tc;
%             diag =  chset(d,h,'diag_1');


            % read molar co2 if it is available

            if ~isempty(chset(d,h,'co2_molar','index'))
                co2 = chset(d,h,'co2_molar');
                h2o = chset(d,h,'h2o_molar');
            elseif ~isempty(chset(d,h,'co2(1)','index'))
                co2 = chset(d,h,'co2(1)')*.14-100;  % convert voltage to concentration
                h2o = chset(d,h,'h2o(1)')*.008-10;  % convert voltage to concentration
            end

            % store time
            [it] = despike(Ts,6,223,350,'Ts'); % changed upper limit from 323 to 350 on June 23 2008, since desert is hotter than canada
            [ic] = despike(co2,6,300,600,'co2');
            %[ih] = despike(h2o,6,0,40,'h2o');
            [ih] = despike(h2o,6,-10,40,'h2o'); % the h2o despike was set too high

            IRGADIAG = (ic & ih);
            %NGOOD_IRGA=length(find(IRGADIAG));

            %fix up the SJER irga h2o
            if iSite == 8;
               [h2o] = SJERH2Oadjust(TIME, h2o);
            end
            
            %fix up the shorthair irga h2o
            if iSite == 9;
               [h2o] = shorthairH2Oadjust(TIME, h2o);
            end
                    
            % Process sonic data if it exists
            % calculate momentum flux
            [UVWROT,SONDIAG,UVWTMEAN_,THETA_,UVWTVAR_,COVUVWT_,UVWTSKEW_,UVWTKURT_,USTAR_,HBUOYANT_] = fluxcsat3(uvwt,diag,'canada');
            NGOOD_SONIC = length(find(SONDIAG));
            NGOOD_SONIC_AND_IRGA = length(find(IRGADIAG&SONDIAG));

            D(176,ii)=NGOOD_SONIC;
            D(177,ii)=NGOOD_SONIC_AND_IRGA;
            
            %Changed so data is directly inserted into designated rows - not appended to the end of a list - awf 5/14/12
            if NGOOD_SONIC > SF*900  %Ensure that there is sufficient data
                if ~isempty(UVWTMEAN_)
                    D(4:7,ii) =  UVWTMEAN_  ; 
                    D(8:11,ii)  =  UVWTVAR_   ; 
                    D(12:15,ii) =  UVWTSKEW_  ; 
                    D(16:19,ii) =  UVWTKURT_  ; 
                    D(21,ii)      =  THETA_     ; 
                    D(22:27,ii)  =  COVUVWT_   ; 
                    D(20,ii)      =  USTAR_     ; 
                    D(30,ii)   =  HBUOYANT_  ; 
                end
            else   % added else clause to deal with situation that SONIC spikes > 900
                D(4:7,ii) =  -999 * ones(4,1)  ; 
                D(8:11,ii)  =    -999 * ones(4,1); 
                D(12:15,ii) =   -999 * ones(4,1) ; 
                D(16:19,ii)  =   -999 * ones(4,1) ; 
                D(21,ii)     =    -999 * ones(1,1); 
                D(22:27,ii)  =    -999 * ones(6,1); 
                D(20,ii)      =    -999 * ones(1,1); 
                D(30,ii)  =   -999 * ones(1,1) ; 
                
                disp(['NGOOD SONIC TOO LOW! = ' NGOOD_SONIC '...filling SONIC REQD DATA with NaNs']);
            end

            % Process licor data if it exists
            if NGOOD_SONIC_AND_IRGA > SF*900   
                
                %removed shaw 5/31/2012 - awf
                
                 %new Tau finder 
                 TauandDelay0 = 0
                %if TauandDelay == 1 && hh <= 2 || TauandDelay == 1 && hh >= 15;%excludes nightime searches for Tau
                %    TauandDelay0 = 1;
                %else
                %    TauandDelay0 = 0;
                %end

                %this is the Tau that is used
                PFlux.tauco2 = tauRefc(1);       TAUC(2,ii) = tauRefc(1); 
                PFlux.tauh2o = tauRefh(1);       TAUH(2,ii) = tauRefh(1); 
                PFlux.delayco2 = delayRefc(1);   DELAYC(2,ii) = delayRefc(1); 
                PFlux.delayh2o = delayRefh(1);   DELAYH(2,ii) = delayRefh(1); 


                %[CO2_,H2O_,RHOM_,FCO2_,FH2O_,HLATENT_,HSENSIBLE_,GAINS_, Ntemp_, Nco2_, Nh2o_]=FluxClosedPath(UVWROT(3,:),Ts,SONDIAG,co2,h2o,IRGADIAG,PFlux,sitePress(iSite));

                
                if TauandDelay == 1
                    %this helps you find the correct Tau and delay
                    [CO2_,H2O_,RHOM_,FCO2_,FH2O_,HLATENT_,HSENSIBLE_,GAINS_, Ntemp_, Nco2_, Nh2o_, Tau_c_max_, Tau_h2o_max_, TauCorrRecord_]=FluxClosedPathDelayandTau3TestsNoShaw(UVWROT(3,:),Ts,SONDIAG,co2,h2o,IRGADIAG,PFlux,sitePress(iSite), TIME, TauandDelay0);
                    pause
                else
                    %once you have found the delay and Tau use this path
                    [CO2_,H2O_,RHOM_,FCO2_,FH2O_,HLATENT_,HSENSIBLE_,GAINS_, Ntemp_, Nco2_, Nh2o_]=FluxClosedPathNoShaw(UVWROT(3,:),Ts,SONDIAG,co2,h2o,IRGADIAG,PFlux,sitePress(iSite));
                    Tau_c_max_=NaN;
                    Tau_h2o_max_=NaN;
                end
                
                %this is the calculated tau
                TAUC(1,ii)   = Tau_c_max_; 
                TAUH(1,ii)   = Tau_h2o_max_; 
                DELAYC(1,ii) = Nco2_/4; %convert into seconds from n shifts
                DELAYH(1,ii) = Nh2o_/4; 
                    
                
                %Insert Taus into D
                D(49:50,ii)=TAUC(:,ii);
                D(51:52,ii)=TAUH(:,ii);
                D(53:54,ii)= DELAYC(:,ii);
                D(55:56,ii)= DELAYH(:,ii);
                
                
                if ~isempty(FCO2_)
                    D(31:37,ii)  =  CO2_  ; 
                    D(45,ii) =  FCO2_  ; 
                    D(38:44,ii)  =  H2O_   ; 
                    D(46,ii) =  FH2O_  ; 
                    D(47:48,ii) =  RHOM_  ; 
                    D(29,ii)  =  HLATENT_    ; 
                    D(28,ii)=  HSENSIBLE_    ; 
                    D(57,ii)   =  GAINS_(1); 
                    D(58,ii)   =  GAINS_(2);
                    %for delay and Tau
                    D(180,ii)   =  Ntemp_;
                    D(181,ii)   =  Nco2_;
                    D(182,ii)   =  Nh2o_;
                    D(183,ii)   =  Tau_c_max_;
                    D(184,ii)   =  Tau_h2o_max_;
                    %end delay and Tau
                    if TauandDelay == 1
                        TCR=[TCR;TauCorrRecord_];
                    end
                end

            else
                D(31:37,ii)  =  -999 * ones(7,1)  ; 
                D(38:44,ii)  = -999 * ones(7,1)  ; 
                D(45,ii) =  -999 * ones(1,1); 
                D(46,ii) =  -999 * ones(1,1)  ; 
                D(47:48,ii) = -999 * ones(2,1)  ; 
                D(29,ii)  = -999 * ones(1,1)  ; 
                D(28,ii)= -999 * ones(1,1)  ; 
                D(57,ii)   = -999 * ones(1,1)  ; 
                D(58,ii)   = -999 * ones(1,1)  ; 
                D(49:50,ii)   = -999 * ones(2,1)  ; 
                D(51:52,ii)   = -999 * ones(2,1)  ; 
                D(53:54,ii) = -999 * ones(2,1)  ; 
                D(55:56,ii) = -999 * ones(2,1)  ; 
                %for delay and Tau
                D(180,ii)   =  NaN;
                D(181,ii)   =  NaN;
                D(182,ii)   =  NaN;
                D(183,ii)   =  NaN;
                D(184,ii)   =  NaN;
                %end delay and Tau
                disp(['NGOOD SONIC&IRGA TOO LOW! = ' NGOOD_SONIC_AND_IRGA ' ...filling SONIC+IRGA REQD DATA with NaNs']);
                %pause;
            end
        else
            %if no fast data is there, then fill with NaNs
            D(4:58, ii) =-999;
            disp(' no data returned from reading file');
            
        end  % executes if tower top is working     %L3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %       Now deal with the slow data
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %The output of the data logger data is not consistent with time.
        %Therefore, we have to search the header and put the data into a
        %consistent format. - awf 5/14/12
        
        %First, join the slow data tables together and make sure we only
        %have 1 interval
        
        d30=[];
        h30=[];

        %(A) data logger data
        if ~isempty(dflux)
            dflux(isnan(dflux))=-999;

            if size(dflux,2) >1
                s=sprintf('\n MORE THAN ONE INTERVAL RETURNED FOR THE FLUX DATA (CHECK FOR INTERVAL [ %8.4f ]) \n',JDInterval(ii)+1/48);
                fprintf(1,'\n %s \n',s);
                %RunInfo=[RunInfo s];
            end

            d30 = dflux(1:end,1);
            h30 = char(hflux(1:end,:));
        end

        % (B) Append the data from the slow variance files
        % (eg 2001_07_27_0030_variance.TOB) files

        if ~isempty(dvariance)
            dvariance(isnan(dvariance))=-999;

            if size(dvariance,2) >1
                s=sprintf('\n MORE THAN ONE INTERVAL RETURNED FOR THE VARIANCE DATA (CHECK FOR INTERVAL [ %8.4f ]) \n',JDInterval(ii)+1/48);
                fprintf(1,'\n %s \n',s);
                %RunInfo=[RunInfo,s];
            end

            d30 = [d30;dvariance(3:end,1)];
            h30 = strvcat(h30,char(hvariance(3:end,:))); %#ok<VCAT>
        end
        
        %(C) Search the datalogger and variance data, find the variable of interest, 
        %and put data into the D output data matrix
        
        if ~isempty(dflux) && ~isempty(dvariance) 
            %record datalogger data
            %SLOW_D30= [SLOW_D30, d30(:,1)];

            %first fill the D interval with -999;
            D(59:179,ii)=-999;

            %put the headers into lower case to facilitate header string
            %agreement
            lower_MVL_Universal = lower(MVL_Universal);
            lower_h30 = lower(h30);

            %LOC is lower_MVL_Universal long with numbers indexing the
            %last row (in our case, the only row) in h30/d30 with the corresponding data
            %and 0 where there is no agreement
            [TF,LOC] = ismember(lower_MVL_Universal,lower_h30);

            good=LOC>0; 
            rowindex=LOC(good); 
            D(TF, ii)= d30(rowindex,1);
        else
            D(59:179,ii)=-999;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% save data once per day
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if abs( HHMM - 1200 ) < 10 || ii==i2Proc(end)
            
            D(isnan(D)) =-999;
            disp(['Saving to ' fileout_one_array ' and to ' fileout_one_arraybak]);
            eval(['save ' fileout_one_array ' D MVL_Universal RunInfo']);
            eval(['save ' fileout_one_arraybak ' D MVL_Universal RunInfo']);
            disp('Finished Saving');
        
        end   % for each half hour interval
    
    end % of processing each site

disp('Time to process data is ');
toc

diary off
end
return