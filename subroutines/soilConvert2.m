function Soil_data = soilConvert(sites2Proc)

var_defs();

global sites soilDir dirSep towerYearStart;

%%
% Output header format:
% Cols  1-10: Towerdate,towerDay,towerHour,towerMinute,record,batt_volt,RefTemp_C,Soil_M(1),Soil_M(2),Soil_M(3),
% Cols 11-20: Soil_M(4),Fuel_M(1),Fuel_M(2),Soil_T(1),Soil_T(2),Soil_T(3),Soil_T(4),Box_T,LWS(1),LWS(2),
% Cols 21-30: LWS(3),StartT_C(1),StartT_C(2),StartT_C(3),StartT_C(4),StartT_C(5),StartT_C(6),T_1s_C(1),T_1s_C(2),T_1s_C(3),
% Cols 31-40: T_1s_C(4),T_1s_C(5),T_1s_C(6),T_30s_C(1),T_30s_C(2),T_30s_C(3),T_30s_C(4),T_30s_C(5),T_30s_C(6),DelT_C(1),
% Cols 41-47: DelT_C(2),DelT_C(3),DelT_C(4),DelT_C(5),DelT_C(6),AirT,SnowD
%%

for iSite=sites2Proc
    siteName = char(sites(iSite));
    
    soilSiteDir = [soilDir siteName dirSep];
    soilFiles = dir([soilSiteDir '*.dat']);
    %Soil_data = NaN * ones(1,46);
    Soil_data = [];
   
    % read in the files
    for ii=1:length(soilFiles)
        soilFile = [soilSiteDir soilFiles(ii).name];
        disp(['Reading ' soilFile]);

        
        %% some data are pure doubles in a csv format, so see if it's one
        % of those first
        try
            fileData = csvread(soilFile);
            hours = floor(fileData(:,4)/100);
            minutes = (fileData(:,4)/100-floor(fileData(:,4)/100))*100;
            soilDates = datenum(fileData(:,2),0,fileData(:,3),hours,minutes,0)-towerYearStart(iSite);
            soilData = [soilDates fileData(:,2:end)];
            soilData = soilData';
            
            Soil_M = NaN * ones(4,size(soilData,2));
            Fuel_M = NaN * ones(2,size(soilData,2));
            Soil_T = NaN * ones(4,size(soilData,2));
            LWS    = NaN * ones(3,size(soilData,2)); 
            StartTemp_C  = NaN * ones(6,size(soilData,2));
            Temp_1sec_C  = NaN * ones(6,size(soilData,2));
            Temp_30sec_C = NaN * ones(6,size(soilData,2));
            DeltaT_C     = NaN * ones(6,size(soilData,2));
            box_T        = NaN * ones(1,size(soilData,2));
            airT         = NaN * ones(1,size(soilData,2));
            snowD        = NaN * ones(1,size(soilData,2));
            
            % towerDates	YYYY  DD HHMM SS batt_volt	RefTemp_C	Soil_M(1)	Soil_M(2)	Soil_M(3)	Soil_M(4)	Soil_T(1)	Soil_T(2)	Soil_T(3)	Soil_T(4)	StartTemp_C(1)	StartTemp_C(2)	StartTemp_C(3)	StartTemp_C(4)	Temp_1sec_C(1)	Temp_1sec_C(2)	Temp_1sec_C(3)	Temp_1sec_C(4)	Temp_30sec_C(1)	Temp_30sec_C(2)	Temp_30sec_C(3)	Temp_30sec_C(4)	DeltaT_C(1)	DeltaT_C(2)	DeltaT_C(3)	DeltaT_C(4)	AirT	SnowD
            towerDate = soilData(1,:);
            towerDay  = floor(towerDate);
            [a b c hours mins d] = datevec(towerDate); %#ok<*ASGLU>
            batt_volt = soilData(6,:);
            RefTemp_C = soilData(7,:);
            Soil_M(1:4,:) = soilData(8:11,:);
            Soil_T(1:4,:) = soilData(12:15,:);
            StartTemp_C(1:4,:)  = soilData(16:19,:);
            Temp_1sec_C(1:4,:)  = soilData(20:23,:);
            Temp_30sec_C(1:4,:) = soilData(24:27,:);
            DeltaT_C(1:4,:)     = soilData(28:31,:);
            airT                = soilData(32,:);
            snowD               = soilData(33,:);
            
        catch err
            if (strcmp(err.identifier,'MATLAB:textscan:handleErrorAndShowInfo'))
                %% read in the data from the *dat files
                fid = fopen(soilFile,'r');
                lineArray = {};
                lineIndex = 1;
                nextLine  = fgetl(fid);
                while ~isequal(nextLine,-1)
                    lineArray{lineIndex} = nextLine;
                    lineIndex = lineIndex+1;
                    nextLine  = fgetl(fid);
                end
                fclose(fid);

                lineArray = lineArray(1:lineIndex-1);
                for iLine = 1:lineIndex-1
                    lineData = textscan(lineArray{1,iLine}, '%s', 'Delimiter', ',');
                    lineData = lineData{1};
                    %if strcmp(lineArray{iLine}(end), ',')
                    %    lineData{end+1} = '';
                    %end
                    lineArray(1:numel(lineData),iLine) = lineData;
                end
                lineDoubles = cellfun(@(s) {str2double(s)},lineArray);

                % find the beginning of the first data row
                recRow = cell2mat(lineDoubles(2,:));
                firstData = find(~isnan(recRow),1);

                % convert the date stamps to tower dates
                soilDates = [];
                for jj = firstData:size(lineArray,2) %for jj = firstData:length(lineArray)
                    soilDates = [soilDates datenum(lineArray(1,jj)) - towerYearStart(iSite)];
                end

                soilData = cell2mat(lineDoubles);
                soilData = soilData(:,firstData:end);
                soilData(1,:) = soilDates;

                %% there are three TOA header types right now, luckily each of them
                %  is a different length so we can tell the difference
                % TIMESTAMP	RECORD	batt_volt	RefTemp_C	Soil_M(1)	Soil_M(2)	Soil_M(3)	Soil_M(4)	Fuel_M(1)	Fuel_M(2)	Soil_T(1)	Soil_T(2)	Soil_T(3)	Soil_T(4)	Box_T	LWS(1)	LWS(2)	LWS(3)																								
                % TMSTAMP	RECNBR	batt_volt	RefTemp_C	Soil_M(1)	Soil_M(2)	Soil_M(3)	Soil_M(4)	Soil_T(1)	Soil_T(2)	Soil_T(3)	Soil_T(4)	StartTemp_C(1)	StartTemp_C(2)	StartTemp_C(3)	StartTemp_C(4)	Temp_1sec_C(1)	Temp_1sec_C(2)	Temp_1sec_C(3)	Temp_1sec_C(4)	Temp_30sec_C(1)	Temp_30sec_C(2)	Temp_30sec_C(3)	Temp_30sec_C(4)	DeltaT_C(1)	DeltaT_C(2)	DeltaT_C(3)	DeltaT_C(4)	AirT	SnowD												
                [leng w] = size(soilData);

                %% output format:
                % towerDate towerDay hours mins batt_volt RefTemp_C Soil_M(1) Soil_M(2) Soil_M(3) Soil_M(4) Fuel_M(1) Fuel_M(2) Soil_T(1) Soil_T(2) Soil_T(3) Soil_T(4) Box_T LWS(1) LWS(2) LWS(3) StartT_C(1) StartT_C(2) StartT_C(3) StartT_C(4) StartT_C(5) StartT_C(6) T_1s_C(1) T_1s_C(2) T_1s_C(3) T_1s_C(4) T_1s_C(5) T_1s_C(6) T_30s_C(1) T_30s_C(2) T_30s_C(3) T_30s_C(4) T_30s_C(5) T_30s_C(6) DelT_C(1) DelT_C(2) DelT_C(3) DelT_C(4) DelT_C(5) DelT_C(6) airT snowD																										
                Soil_M = NaN * ones(4,size(soilData,2));
                Fuel_M = NaN * ones(2,size(soilData,2));
                Soil_T = NaN * ones(4,size(soilData,2));
                LWS    = NaN * ones(3,size(soilData,2)); 
                StartTemp_C  = NaN * ones(6,size(soilData,2));
                Temp_1sec_C  = NaN * ones(6,size(soilData,2));
                Temp_30sec_C = NaN * ones(6,size(soilData,2));
                DeltaT_C     = NaN * ones(6,size(soilData,2));
                box_T        = NaN * ones(1,size(soilData,2));
                airT         = NaN * ones(1,size(soilData,2));
                snowD        = NaN * ones(1,size(soilData,2));


                if leng == 18
                % TIMESTAMP	RECORD	batt_volt	RefTemp_C Soil_M(1)	Soil_M(2) Soil_M(3)	Soil_M(4)	Fuel_M(1)	Fuel_M(2)	Soil_T(1)	Soil_T(2)	Soil_T(3)	Soil_T(4)	Box_T	LWS(1)	LWS(2)	LWS(3)
                    towerDate = soilData(1,:);
                    towerDay  = floor(towerDate);
                    [a b c hours mins d] = datevec(towerDate);
                    batt_volt = soilData(3,:);
                    RefTemp_C = soilData(4,:);
                    Soil_M(1:4,:) = soilData(5:8,:);
                    Fuel_M(1:2,:) = soilData(9:10,:);
                    Soil_T(1:4,:) = soilData(11:14,:);
                    box_T         = soilData(15,:);
                    LWS(1:3,:)    = soilData(16:18,:);
                elseif leng == 30
                % TMSTAMP	RECNBR	batt_volt	RefTemp_C	Soil_M(1)	Soil_M(2)	Soil_M(3)	Soil_M(4)	Soil_T(1)	Soil_T(2)	Soil_T(3)	Soil_T(4)	StartTemp_C(1)	StartTemp_C(2)	StartTemp_C(3)	StartTemp_C(4)	Temp_1sec_C(1)	Temp_1sec_C(2)	Temp_1sec_C(3)	Temp_1sec_C(4)	Temp_30sec_C(1)	Temp_30sec_C(2)	Temp_30sec_C(3)	Temp_30sec_C(4)	DeltaT_C(1)	DeltaT_C(2)	DeltaT_C(3)	DeltaT_C(4)	AirT	SnowD
                    towerDate = soilData(1,:);
                    towerDay  = floor(towerDate);
                    [a b c hours mins d] = datevec(towerDate);
                    batt_volt = soilData(3,:);
                    RefTemp_C = soilData(4,:);
                    Soil_M(1:4,:) = soilData(5:8,:);
                    Soil_T(1:4,:) = soilData(9:12,:);
                    StartTemp_C(1:4,:) = soilData(13:16,:);
                    Temp_1sec_C(1:4,:) = soilData(17:20,:);
                    Temp_30sec_C(1:4,:) = soilData(21:24,:);
                    DeltaT_C(1:4,:)     = soilData(25:28,:);
                    airT                = soilData(29,:);
                    snowD               = soilData(30,:);
                elseif leng == 33

                elseif leng == 42
                % TIMESTAMP	RECORD	batt_volt	RefTemp_C	Soil_M(1)	Soil_M(2)	Soil_M(3)	Soil_M(4)	Fuel_M(1)	Fuel_M(2)	Soil_T(1)	Soil_T(2)	Soil_T(3)	Soil_T(4)	Box_T	LWS(1)	LWS(2)	LWS(3)	StartTemp_C(1)	StartTemp_C(2)	StartTemp_C(3)	StartTemp_C(4)	StartTemp_C(5)	StartTemp_C(6)	Temp_1sec_C(1)	Temp_1sec_C(2)	Temp_1sec_C(3)	Temp_1sec_C(4)	Temp_1sec_C(5)	Temp_1sec_C(6)	Temp_30sec_C(1)	Temp_30sec_C(2)	Temp_30sec_C(3)	Temp_30sec_C(4)	Temp_30sec_C(5)	Temp_30sec_C(6)	DeltaT_C(1)	DeltaT_C(2)	DeltaT_C(3)	DeltaT_C(4)	DeltaT_C(5)	DeltaT_C(6)
                    towerDate = soilData(1,:);
                    towerDay  = floor(towerDate);
                    [a b c hours mins d] = datevec(towerDate);
                    batt_volt = soilData(3,:);
                    RefTemp_C = soilData(4,:);
                    Soil_M(1:4,:) = soilData(5:8,:);
                    Fuel_M(1:2,:) = soilData(9:10,:);
                    Soil_T(1:4,:) = soilData(11:14,:);
                    box_T         = soilData(15,:);
                    LWS(1:3,:)    = soilData(16:18,:);
                    StartTemp_C(1:6,:) = soilData(19:24,:);
                    Temp_1sec_C(1:6,:) = soilData(25:30,:);
                    Temp_30sec_C(1:6,:) = soilData(31:36,:);
                    DeltaT_C(1:6,:)     = soilData(37:42,:);
                else
                    error('*** Error: File type not identified. Reading soil data aborted.');
                end
            else
                rethrow(err);
            end
        end
        %% aggregate data into one big array
        Soil_data = [towerDate' towerDay' hours' mins' batt_volt' RefTemp_C' Soil_M' Fuel_M' Soil_T' box_T' LWS' StartTemp_C' Temp_1sec_C' Temp_30sec_C' DeltaT_C' airT' snowD'; Soil_data];
    end
    %% delete duplicate entries, put in placeholders for missing times
    Soil_data = sortrows(Soil_data);
    Soil_data(Soil_data == -7999) = NaN;
    Soil_data(Soil_data(:,1)<0,:) = [];  % delete weird rows with negative timestamps
    [junk,jj,kk] = unique(Soil_data(:,1),'first');
    Soil_data = Soil_data(jj,:);
    allTimes = Soil_data(1,1):1/48:Soil_data(end,1);
    soilFill = NaN * ones(length(allTimes),size(Soil_data,2));
    soilFill(:,1) = allTimes;
    [junk, ia, ib] = intersect(soilFill(:,1),Soil_data(:,1));
    soilFill(ia,:) = Soil_data(ib,:);    
 
    C = setdiff(Soil_data(1,5):1:Soil_data(end,5),Soil_data(:,5));
    if length(C) > 1
        disp(['*** Warning: missing ' num2str(length(C)) ' data points in dataset.']);
    else
        disp(['Dataset complete from day ' num2str(Soil_data(1,1)) ' to ' num2str(Soil_data(end,1))]);
    end
    outputFile = [soilSiteDir siteName '_soil.csv'];
    dlmwrite(outputFile,Soil_data,'precision',10);
end

end