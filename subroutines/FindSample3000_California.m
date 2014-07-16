function [FileName,ISAMPLE]=FindSample3000_California(Directory,Time2Find,MaxSecOff)

% looks in cr3000 output files to find the sample nearest to
% Time2Find, but not more than MaxSecOff away from it. Time2Find
% is in days since Jan 1, 2006

% going to have to deal with cases where two files were created
% for the same day, meaning the program or something changed.

% make sure Directory ends in a backslash


if isunix
    DirectorySeparator='/';
else
    DirectorySeparator='\';
end


if ~strcmp(Directory(end),DirectorySeparator);Directory=[Directory DirectorySeparator];end

% put time 2 find into seconds to find

 Seconds2Find =  (Time2Find+(datenum('00-jan-2008')-datenum('01-jan-1990')))*24*3600;   %% amsm 10/18/04
%Seconds2Find =  (Time2Find+(datenum(2001,1,0)-datenum('01-jan-1990')))*24*3600;
% put file names into Seconds

%  get a list of the data files in the directory

files = dir([Directory]);
%put in some logic to avoid unix measuring the width of the files array
%incorrectly
if ~isempty(findstr('ts_data',Directory)),
    fnamelength=37;
    
    files = char(files.name);
    files=files(:,1:fnamelength);
    files_lower_names = lower(files);%this variable created for .tob comparison

    itablefile=strmatch('tob1',files_lower_names(:,1:4));
    clear files_lower_names;
    files = files(itablefile,:);
    % disp(files); pause;
    % sort the files by date

    FileYear  = str2num(files(:,[ 19:22]));
    FileMonth = str2num(files(:,[ 24  25]));
    FileDay   = str2num(files(:,[ 27 28]));
    FileHour  = str2num(files(:,[30 31]));
    FileMin   = str2num(files(:,[32 33]));
    
    filedates = str2num(files(:,[19:22 24 25 27 28 30:33]));
elseif ~isempty(findstr('flux',Directory)),
    fnamelength=34;
    
        files = char(files.name);
    files=files(:,1:fnamelength);
    files_lower_names = lower(files);%this variable created for .tob comparison

    itablefile=strmatch('tob1',files_lower_names(:,1:4));
    clear files_lower_names;
    files = files(itablefile,:);
    % disp(files); pause;
    % sort the files by date

    FileYear  = str2num(files(:,[ 16:19]));
    FileMonth = str2num(files(:,[ 21  22]));
    FileDay   = str2num(files(:,[ 24 25]));
    FileHour  = str2num(files(:,[27 28]));
    FileMin   = str2num(files(:,[29 30]));
    
    filedates = str2num(files(:,[16:19 21 22 24 25 27:30]));
elseif ~isempty(findstr('variance',Directory)),
    fnamelength=38;
    
    files = char(files.name);
    files=files(:,1:fnamelength);
    files_lower_names = lower(files);%this variable created for .tob comparison

    itablefile=strmatch('tob1',files_lower_names(:,1:4));
    clear files_lower_names;
    files = files(itablefile,:);

    % sort the files by date

    FileYear  = str2num(files(:,[ 20:23]));
    FileMonth = str2num(files(:,[ 25  26]));
    FileDay   = str2num(files(:,[ 28 29]));
    FileHour  = str2num(files(:,[31 32]));
    FileMin   = str2num(files(:,[33 34]));

    
    filedates = str2num(files(:,[20:23 25 26 28 29 31:34]));
else
    disp('Directory does not match file type in Findsample')
end

%now reduce the width of the files array to the length of the expected
%filename


%FileTimes = datenum(FileYear,FileMonth,FileDay,FileHour,FileMin,0)-datenum('01-jan-2008');
FileTimes = datenum(FileYear,FileMonth,FileDay,FileHour,FileMin,0)-datenum('00-jan-2008');

% sort the files
[FileTimes,index] = sort(FileTimes);



% find the last file that began before the desired interval   
%%%% PROBLEM FOUND HERE 11/03/05. Time2Find for a whole number (eg 1643)
%%%% ifile returns an index that points to a file that contains the data
%%%% for the previous day. The reason why is that there is a precision
%%%% issue. ifile for 1643.0000 should be 1436 NOT 1435 but

%%%%    Eg, for 1643.00000, Time2Find is calcd as 1643.020833333333
%%%%             but FileTimes(1436)  is calcd as 1643.020833333372
%%%% so it is more than Time2Find. We need to round off to 6 decimal places

FileTimes_round = 1e-6 * fix(FileTimes*1e6);
Time2Find_round = 1e-6 * fix(Time2Find*1e6);

ifile = max( find( FileTimes_round <= Time2Find_round) );


if ~isempty(ifile)
    %disp('about to invoke read3000 on line 60 of Find Sample');
    %disp([Directory files(ifile,:)]);
    %disp('paused at line 62 of findsampl before going into read3000');%pause;
    [h,d]=read3000([Directory files(ifile,:)],char({'SECONDS';'NANOSECONDS'}),[],0);

    Seconds    = d(1,:) + d(2,:)/1e9 ;
    ISAMPLE    = find( abs(Seconds2Find-Seconds) == min(abs(Seconds2Find-Seconds)) );
    SecondsOff = abs(Seconds2Find-Seconds(ISAMPLE));

    if SecondsOff > MaxSecOff
        FileName = [];
        ISAMPLE = [];
    else
        FileName = files(ifile,:);
    end

else

    FileName = [];
    ISAMPLE = [];

end


return




% convert the desired timestamp to find into similar format
% as the file naming convention

[Y,M,D,H,MI,S] = datevec(Time2Find+datenum('01-jan-2008'));

if M <10; Month=['0' int2str(M)]; else Month=int2str(M);end
if D <10; Day  =['0' int2str(D)]; else Day  =int2str(D);end
if H <10; Hour =['0' int2str(H)]; else Hour =int2str(H);end
if MI<10; Min  =['0' int2str(MI)];else Min  =int2str(MI);end

% convert 0 hours to 2400 hours

str2num([Hour Min])

%if str2num([Hour Min])==0
%  Hour='24';Min='00';
%end


filestr=[int2str(Y) '_' Month '_' Day]


%disp('paused');%pause;

    [filedates,index] = sort(filedates);

    files = files(index,:);

    ich=strmatch(filestr,files)


if ~isempty(ich)

    % if there were more than 1 files found, check the minutes to see
    % which one should be opened

    if length(ich)>1

        % get hours and minutes for all files that have the same
        % starting date

        HHMM = str2num(files(ich,12:15))

        ifile = max( find( str2num([Hour Min]) >= HHMM ) )

        if ~isempty(ifile)
            ich=ich(ifile)
        end

    end

    % search the file for the nearest sample to
    % Time2Find. In SECONDS, Time2Find is given by our Canada
    % dataset reference start-time of Jan-01-2006, plus the
    % number of seconds accrued from 1990 to 2006

    Seconds2Find =  (Time2Find+(datenum('01-jan-2008')-datenum('01-jan-1990')))*24*3600;

    % consider the case where there are more than 1 file in ich (so
    % ifile didn't return a clear result).  open both of the files
    % and find the one that has the value closest to the desired
    % value

    for jj=1:length(ich)

        % read in the time vector for this file

        [h,d]=read3000([Directory files(ich(jj),:)],char({'SECONDS';'NANOSECONDS'}),[],0);

        Seconds        = d(1,:) + d(2,:)/1e9;
        isample(jj)    = find( abs(Seconds2Find-Seconds) == min(abs(Seconds2Find-Seconds)) )
        SecondsOff(jj) = abs(Seconds2Find-Seconds(isample(jj)))
        SecondsOff2(jj) = Seconds2Find-Seconds(isample(jj))

        rem(Seconds(1)/24/3600,1)

    end

    icloser = min(SecondsOff)

    %
    %   ISAMPLE       = find( abs(Seconds2Find-Seconds) == min(abs(Seconds2Find-Seconds)) )


    if abs(Seconds2Find-Seconds(isample(icloser))) > MaxSecOff
        Filename = [];
        ISAMPLE = [];
    else
        FileName = files(ich(icloser),:)
        ISAMPLE = isample(icloser)
    end

    %disp('paused');%pause;


else

    FileName = [];
    ISAMPLE = [];

end


return


