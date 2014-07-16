%function mergeV2(sites2Proc)
sites2Proc=6;
%LOOK at the significant figures on the xlsread for the goes data...merge
%shows fewer digits and this may explain the historic time offsets - awf

% current version 9/23/10 aek
path(path, 'C:\towerData\ProcessingScripts\subroutines');

%code to get this to run through
goesRootDir='C:\towerData\goes\';
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           ~~~~  Part 1. Initialilizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global sites iSite VIsites towerYearStart MVL_Universal
% global siteAlt
global mergedRootDir inputRootDir fastRootDir
% flux processing parameters
global procInt
%%
var_defs();  
Day = date;

%removed a note about old and new headers.  We are going to make headers consistent. - awf 5/2/2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    2. File Management Section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Deleted importing file headers from excel - awf 5/2/2012

diary_filename = [mergedRootDir 'merge_log_' Day];
%diary(diary_filename);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    3. Establish Headers and Labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3A    Read GOES Headers either from matlab file or from spreadsheet
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%   Start Merging GOES DATA with PROCESSED (ONE ARRAY DATA) %%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');


for iSite=sites2Proc
    % 3B    Read non-old GOES HEaders either from matlab file or from spreadsheet -
    % need to deal with two separate headers for Sierra and non-Sierra data
            
    %rework data files in GOES.excel so all datasets are structured the same
    %==>removed novel headers - awf 5/2/2012
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    5. Start Site Loop for Merging
    %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    siteName = char(sites(iSite));
    if exist(mergedRootDir, 'dir') ~= 7
        mkdir(mergedRootDir);
    end

    %fout = [mergedRootDir siteName '_MRG'];
    fout = [mergedRootDir siteName '_MRG'];
    fouthold= [mergedRootDir 'Hold' siteName '_MRG'];
    
    disp('----------------------------------------------------------------');
    disp(['Site: ' siteName]);
    disp('----------------------------------------------------------------');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   5.1 Manage GOES Data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% rewrite this - aek
    % The Goes Spreadsheets are periodically copied from Gregs computer (manually) and the Data pages are converted to
    % comma-delimited text files (*.csv). In this section we read the
    % textfiles (or after an initial read, read the matlab version of these
    % files) add new time series, and prepare the GOES data for merging with
    % the Processed Data (ie the output from FSC_calfiornia, saved in one array
    % files)

    %Excel GOES.csv were manipulated to be 95 long and headers and data match - awf
    
    % 5.1 A    Read the GOES Data either from the matlab files or from the csv
    % files
    
    %Now can only read data from .csv files - remove mat - awf 5/2/2012
    fin_goes = [goesRootDir siteName '.csv'];

    if exist(fin_goes, 'file') 
        
        [hgoes,dgoes]=goesread(fin_goes);
        
        hgoes = strcat({'GOES_'},hgoes);
        yy = dgoes(2,:);
        jday = dgoes(3,:);
        HH = dgoes(4,:);
        MM = dgoes(5,:);


        MLDT_GOES = datenum(yy,1,0) + jday + HH/24 + MM/1440; % to bump back to tower time
        EXPDAY_GOES = MLDT_GOES - towerYearStart(iSite);
    
    else
        %no goes ==> stop
        disp('missing GOES')
        pause
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    5.2 Manage Processed Data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fin_proc = [fastRootDir siteName '\proc\Fluxes_' siteName '_one_array'];
    eval(['load ' fin_proc ';']);

    %Time for the processed data
    EXPDAY_PROC = D(2,:);
    
    %MLDT_PROC = EXPDAY_PROC + towerYearStart(iSite);

    %We only want the first 180 rows of D.  The other rows contained
    %data that was not needed or duplicate data
    D=D(1:179,:);
    MVL_Universal=cellstr(MVL_Universal(1:179,1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    5.3 Create New Array of Times that incorportes both GOES and PROC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    First_Valid_GOES_Day = min(EXPDAY_GOES(isfinite(EXPDAY_GOES) & EXPDAY_GOES>0));
    Last_Valid_GOES_Day = max(EXPDAY_GOES(isfinite(EXPDAY_GOES) & EXPDAY_GOES>0));

    First_Valid_Proc_Day = min(EXPDAY_PROC(isfinite(EXPDAY_PROC) & EXPDAY_PROC>0));
    Last_Valid_Proc_Day = max(EXPDAY_PROC(isfinite(EXPDAY_PROC) & EXPDAY_PROC>0));

    First_Valid_Day_Either = min([ First_Valid_GOES_Day First_Valid_Proc_Day]);
    Last_Valid_Day_Either = max([Last_Valid_GOES_Day Last_Valid_Proc_Day]);

    length_valid_proc = length(find(isfinite(EXPDAY_PROC) & EXPDAY_PROC>0)); %#ok<NASGU>
    length_valid_goes = length(find(isfinite(EXPDAY_GOES) & EXPDAY_GOES>0)); %#ok<NASGU>

    EXPDAY_MERGE = First_Valid_Day_Either: procInt : Last_Valid_Day_Either ;
    NDATA_Merge = length(EXPDAY_MERGE);

    % now truncate the decimal places to 6 digits
    EXPDAY_MERGE = fix(EXPDAY_MERGE*1e6)*1e-6;
    EXPDAY_PROC = fix(EXPDAY_PROC*1e6)*1e-6;
    EXPDAY_GOES = fix(EXPDAY_GOES*1e6)*1e-6;

    [cix_proc,ix_proc,iy_proc] = intersect(EXPDAY_MERGE,EXPDAY_PROC);
    [cix_goes,ix_goes,iy_goes] = intersect(EXPDAY_MERGE,EXPDAY_GOES);

    disp(['Data ranges from ' datestr(First_Valid_Day_Either + towerYearStart(iSite)) ' to ' ...
        datestr(Last_Valid_Day_Either + towerYearStart(iSite)) ]);
    disp('----------------------------------------------------');
    disp(['Length of MERGE is ' num2str(NDATA_Merge)]);
    disp(['Proccessed Data takes up ' num2str( 100 * length(cix_proc)/NDATA_Merge,3) '% of Merge']);
    disp(['GOES       Data takes up ' num2str( 100 * length(cix_goes)/NDATA_Merge,3) '% of Merge']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    5.4 Define Dimensions and Insert Data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 5.4.1  Work out the positions in DMERGE where the Proc and GOES data go
    NVARS_GOES = size(dgoes,1);
    NVARS_PROC = size(D,1);
    NVARS_MERGE = 1 + NVARS_GOES + NVARS_PROC ;
    PROC_IX_STT = 2;
    PROC_IX_END = PROC_IX_STT + NVARS_PROC -1;
    GOES_IX_STT = PROC_IX_END + 1;
    GOES_IX_END = GOES_IX_STT + NVARS_GOES -1;

    DMERGE = NaN * ones(NVARS_MERGE,NDATA_Merge);

    %5.4.2  Add in Universal Time in EXPDAYS
    DMERGE(1,:) = EXPDAY_MERGE;

    %5.4.3  Insert Data from Processed Data using common timestamps
    DMERGE(PROC_IX_STT:PROC_IX_END , ix_proc ) = D(:,iy_proc);

    %5.4.4  Insert Data from GOES Data using common timestamps
    DMERGE(GOES_IX_STT:GOES_IX_END , ix_goes ) = dgoes(:,iy_goes); %#ok<NASGU>

    %5.4.5  Make a Header File
    HMERGE = [{'EXPDAY_MERGE'} ; MVL_Universal ; cellstr(hgoes)]; %#ok<NASGU>

    %for memory
    clear D dgoes

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    5.5 Add Soil Data - awf
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Load in the soil data from excel - this takes a while
    siteName = char(sites(iSite));
    SoilDir='C:\towerData\SoilMoisture\';
    inputfile = [SoilDir siteName '_Soil.xlsx'];
    [a, H_Soil]=xlsread(inputfile ,'SoilMoistureData', 'A1:AT1'); %load in the header
    H_Soil=H_Soil';
    [a, stime]=xlsread(inputfile ,'SoilMoistureData', 'A:A'); %text time stamp

    Soil_textTime=stime(2:end,1);
    length_D=length(Soil_textTime(:,1))+1;% +1 will read in all the data in the excel sheet
    strlength_D=num2str(length_D);
    
    block1='B2:AT30000';
    block2=['B30001:AT' strlength_D]; 
    %read in data by parts
    [a]=xlsread(inputfile ,'SoilMoistureData', block1);
    Soil_data=a;
    clear a
    [a]=xlsread(inputfile ,'SoilMoistureData', block2);
    Soil_data=[Soil_data; a];
    clear a
    %---------------------------------------------------------

    %reshape the soil data so it looks like merge
    Soil_data=Soil_data';
  
    
    %Get the time from the time stamp -- Omit First row...it was the header
    TimeSoil = datenum(Soil_textTime)-towerYearStart(iSite); 
    TimeSoil=TimeSoil';
    TimeSoil = fix(TimeSoil*1e6)*1e-6;%get the soil time lined up
    
    %Add the new calculated time
    Soil_data=[TimeSoil; Soil_data];
    H_Soil(1,1) = {'TimeSoil_mat_calc'}; %converted TS to matlab calculated time
    
    %an array of NaNs that has n rows = length of the soil header and n columns = width of merge
    %build the array to dump in the numbers
    length_S = length(Soil_data(:,1));
    MergeWidth=size(DMERGE, 2);%This is the width of merge
    
    SoilAdd=ones(length_S, MergeWidth)*NaN;
    
    %Where do the time stamps on the 2 datasets overlap?
    [cix_proc,ix_soil,iy_soil] = intersect(EXPDAY_MERGE, TimeSoil);%data logger
    SoilAdd(:, ix_soil) = Soil_data(:,iy_soil);
     
    %add the soil to merge
    DMERGE=[DMERGE;SoilAdd];
    HMERGE = [HMERGE; H_Soil];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %5.4.6    Save the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    disp(['save ' fout ' DMERGE HMERGE ;']);
    eval(['save ' fout ' DMERGE HMERGE ;']);
    clear hgoes
    
end   % of: for iSites=sites2Proc,
diary off;




