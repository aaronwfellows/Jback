function [FileName,ISAMPLE]=FindSample5000(Directory,Time2Find,MaxSecOff)

global towerYearStart iSite

% looks in cr5000 output files to find the sample nearest to
% Time2Find, but not more than MaxSecOff away from it. Time2Find
% is in days since Jan 1, 2006

% going to have to deal with cases where two files were created
% for the same day, meaning the program or something changed.

% put time 2 find into seconds to find
Seconds2Find = (Time2Find + (towerYearStart(iSite)-datenum(1990,1,1)))*24*3600;

%  get a list of the data files in the directory
files = dir(Directory);

% put in some logic to avoid unix measuring the width of the files array
% incorrectly
if findstr('ts_data',Directory)
    fnamelength=27;
elseif findstr('flux',Directory)
    fnamelength=24;
elseif findstr('variance',Directory)
    fnamelength=28;
else
    disp(Directory);
    disp('FindSample5000: Directory does not match file type in Findsample');
    pause;
end

% now reduce the width of the files array to the length of the expected
% filename
files = char(files.name);
files=files(:,1:fnamelength);
files_lower_names = lower(files);%this variable created for .tob comparison

if isempty(files)
    disp('files variable in FindSample5000.m is empty.  This may cause problems.');
    disp('line 53');
    pause;
end

%whos files;
%disp (Directory)
%pause;
itablefile=strmatch('.tob',files_lower_names(:,end-3:end));
clear files_lower_names;
files = files(itablefile,:);
if isempty(files)
    disp('files variable in FindSample5000.m is empty.  This may cause problems.');
    disp('line 67');
    pause;
end

% sort the files by date
FileYear  = str2num(files(:,1:4));
FileMonth = str2num(files(:,[6  7]));
FileDay   = str2num(files(:,[9 10]));
FileHour  = str2num(files(:,[12 13]));
FileMin   = str2num(files(:,[14 15]));
if (isempty(FileYear) || isempty(FileMonth) || isempty(FileDay) || isempty(FileHour) || isempty(FileMin))
    disp('Time is not going well.  Check FileSample5000.m  Lines 69-73');
    disp(files(1,:));
    disp(files);
    %disp([FileYear FileMonth FileDay FileHour FileMin]);
    %disp([FileYear(1) FileMonth(1) FileDay(1) FileHour(1) FileMin(1)]);
    pause;
end

% if ~isempty(findstr('Sierra_1',Directory))
%     FileTimes = datenum(FileYear,FileMonth,FileDay,FileHour,FileMin,0)-datenum('01-jan-2008');
%  if isempty(FileTimes)
%      disp('FileTimes variable in FindSample5000_California.m is empty.  This may cause problems.');
%  end
% else
%   FileTimes = datenum(FileYear,FileMonth,FileDay,FileHour,FileMin,0)-datenum('01-jan-2006');
%    if isempty(FileTimes)
%      disp('FileTimes variable in FindSample5000_California.m is empty.  This may cause problems.');
%    end
% end

FileTimes = datenum(FileYear,FileMonth,FileDay,FileHour,FileMin,0) - towerYearStart(iSite);

% sort the files
[FileTimes,index] = sort(FileTimes);

files = files(index,:);

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

ifile = find(FileTimes_round <= Time2Find_round, 1, 'last' );

if ~isempty(ifile)
    %disp('about to invoke read5000 on line 60 of Find Sample');
    %disp([Directory files(ifile,:)]);
    %disp('paused at line 62 of findsampl before going into read5000');%pause;
    [Wonky,d]=read5000([Directory files(ifile,:)],char({'SECONDS';'NANOSECONDS'}),[],0);
    %%% files(ifile,:)
    % disp('This Findsample5000 program is starting with a different file from the original _cal progs.');
    % pause;
    Seconds    = d(1,:) + d(2,:)/1e9;
    ISAMPLE    = find( abs(Seconds2Find-Seconds) == min(abs(Seconds2Find-Seconds)) );
    SecondsOff = abs(Seconds2Find-Seconds(ISAMPLE));
    if SecondsOff > MaxSecOff
        FileName = [];
        ISAMPLE = [];
        disp([SecondsOff MaxSecOff]);
        disp(['Seconds off greater than MaxSecOff for file ' files(ifile,:)]);
    else
        FileName = files(ifile,:);
    end

else
    FileName = [];
    ISAMPLE = [];
    disp(['Not finding Sample in file ' files(ifile,:) '.  Check FindSample5000.m']);
end

return
%% return 2/24/12 aek

% %% convert the desired timestamp to find into similar format
% %% as the file naming convention
% [Y,M,D,H,MI,S] = datevec(Time2Find + towerYearStart(iSite) + 1);
%  
% if M <10; Month=['0' int2str(M)]; else Month=int2str(M);end
% if D <10; Day  =['0' int2str(D)]; else Day  =int2str(D);end
% if H <10; Hour =['0' int2str(H)]; else Hour =int2str(H);end
% if MI<10; Min  =['0' int2str(MI)];else Min  =int2str(MI);end
%  
% % % convert 0 hours to 2400 hours
%  
% str2num([Hour Min]);
% 
% if str2num([Hour Min])==0
%   Hour='24';Min='00';
% end
%  
% 
% filestr=[int2str(Y) '_' Month '_' Day];
% 
% 
% %disp('paused');%pause;
% 
% filedates = str2num(files(:,[1:4 6 7 9 10 12:15]));
% 
% [filedates,index] = sort(filedates);
% 
% files = files(index,:);
% 
% ich = strmatch(filestr,files);
% 
% if ~isempty(ich)
% 
%     % if there were more than 1 files found, check the minutes to see
%     % which one should be opened
% 
%     if length(ich)>1
% 
%         % get hours and minutes for all files that have the same
%         % starting date
% 
%         HHMM = str2num(files(ich,12:15));
% 
%         ifile = max( find( str2num([Hour Min]) >= HHMM ) );
% 
%         if ~isempty(ifile)
%             ich=ich(ifile);
%         end
% 
%     end
% 
%     % search the file for the nearest sample to
%     % Time2Find. In SECONDS, Time2Find is given by our Canada
%     % dataset reference start-time of Jan-01-2006, plus the
%     % number of seconds accrued from 1990 to 2006
% 
%     %Seconds2Find =  (Time2Find+(datenum('01-jan-2006')-datenum('01-jan-1990')))*24*3600;
%     Seconds2Find = (Time2Find + (towerYearStart(iSite)-datenum(1990,1,0)))*24*3600;
% 
%     % consider the case where there are more than 1 file in ich (so
%     % ifile didn't return a clear result).  open both of the files
%     % and find the one that has the value closest to the desired
%     % value
% 
%     for jj=1:length(ich)
% 
%         % read in the time vector for this file
% 
%         [h,d]=read5000([Directory files(ich(jj),:)],char({'SECONDS';'NANOSECONDS'}),[],0);
% 
%         Seconds        = d(1,:) + d(2,:)/1e9;
%         isample(jj)    = find( abs(Seconds2Find-Seconds) == min(abs(Seconds2Find-Seconds)) );
%         SecondsOff(jj) = abs(Seconds2Find-Seconds(isample(jj)));
%         SecondsOff2(jj) = Seconds2Find-Seconds(isample(jj));
%         rem(Seconds(1)/24/3600,1);
% 
%     end
% 
%     icloser = find(isample == min(SecondsOff));
% 
%     isample       = find( abs(Seconds2Find-Seconds) == min(abs(Seconds2Find-Seconds)) );
% 
%     if abs(Seconds2Find-Seconds(isample(icloser))) > MaxSecOff
%         Filename = [];
%         ISAMPLE = [];
%     else
%         FileName = files(ich(icloser),:);
%         ISAMPLE = isample(icloser);
%     end
%     %disp('paused');%pause;
% 
% else
%     
%     disp(['Could not find channel ' ich ' in FindSample5000']);
%     FileName = [];
%     ISAMPLE = [];
% 
% end
% return
% 
% 
