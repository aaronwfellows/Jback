function [h,d] = concat3000_california(Directory,channels,interval,N2Skip)

%function [h,d] = concat5000(Directory,channels,interval,N2Skip)
%
% $$$ clear
% $$$ Directory = 'J:\Canada\1930\ts_data\' ;
% $$$ interval =  [205.0058;215.2];
% $$$ channels = [];
% $$$ N2Skip=100;
% $$$ MaxSecOff = 10;

if isunix
    DirectorySeparator='/';
else
    DirectorySeparator='\';
end

if nargin < 4
    N2Skip=0;
end


MaxSecOff = interval(3);
%disp('invoking FindSample30000 from line 35 of concat');

%[FileBeg,SampleBeg]=FindSample3000_california(Directory,interval(1),MaxSecOff);
%[FileEnd,SampleEnd]=FindSample3000_california(Directory,interval(2),MaxSecOff);
[FileBeg,SampleBeg]=FindSample3000_California(Directory,interval(1),MaxSecOff);
[FileEnd,SampleEnd]=FindSample3000_California(Directory,interval(2),MaxSecOff);


if ~isempty(FileBeg) & ~isempty(SampleBeg) & ~isempty(FileEnd) & ~isempty(SampleEnd)

    %  get a list of the data files in the directory


    %put in some logic to avoid unix measuring the width of the files array
    %incorrectly
    if ~isempty(findstr('ts_data',Directory)),
        fnamelength=37;
    files = dir([Directory]);
    files = char(files.name);
    files=files(:,1:fnamelength);
    files_lower_names = lower(files);%this variable created for .tob comparison


    itablefile=strmatch('tob1',files_lower_names(:,1:4));

    files = files(itablefile,:);

    % sort the files by date

    filedates = str2num(files(:,[19:22 24 25 27 28 30:33]));
    [filedates,index] = sort(filedates);

    files = files(index,:);

    ifilebeg=strmatch(FileBeg,files);
    ifileend=strmatch(FileEnd,files);

    d=[];
    h=[];

    elseif ~isempty(findstr('flux',Directory)),
        fnamelength=34;
        
            files = dir([Directory]);
    files = char(files.name);
    files=files(:,1:fnamelength);
    files_lower_names = lower(files);%this variable created for .tob comparison


    itablefile=strmatch('tob1',files_lower_names(:,1:4));

    files = files(itablefile,:);

    % sort the files by date

    filedates = str2num(files(:,[16:19 21 22 24 25 27:30]));
    [filedates,index] = sort(filedates);

    files = files(index,:);

    ifilebeg=strmatch(FileBeg,files);
    ifileend=strmatch(FileEnd,files);

    d=[];
    h=[];
    elseif ~isempty(findstr('variance',Directory)),
        fnamelength=38;
        
            files = dir([Directory]);
    files = char(files.name);
    files=files(:,1:fnamelength);
    files_lower_names = lower(files);%this variable created for .tob comparison


    itablefile=strmatch('tob1',files_lower_names(:,1:4));

    files = files(itablefile,:);

    % sort the files by date

    filedates = str2num(files(:,[20:23 25 26 28 29 31:34]));
    [filedates,index] = sort(filedates);

    files = files(index,:);

    ifilebeg=strmatch(FileBeg,files);
    ifileend=strmatch(FileEnd,files);

    d=[];
    h=[];
    else
        disp('Directory does not match file type in Findsample')
    end

    %now reduce the width of the files array to the length of the expected
    %filename

    if ifilebeg==ifileend
        %disp('ifilebeg=ifileend');
        [h,d] = read3000([Directory DirectorySeparator files(ifilebeg,:)],channels,[SampleBeg SampleEnd],N2Skip);

    else
        %disp('ifilebeg<>ifileend');
        for iday=ifilebeg:ifileend
            %disp('now reading the complete required data set, we are at l62 of concat5000 about to gointo read3000 ');

            if iday==ifilebeg
                [htmp,dtmp] = read3000([Directory DirectorySeparator files(iday,:)],channels,[SampleBeg NaN],N2Skip);
            elseif iday==ifileend
                [htmp,dtmp] = read3000([Directory DirectorySeparator files(iday,:)],channels,[1 SampleEnd],N2Skip);
            else
                [htmp,dtmp] = read3000([Directory DirectorySeparator files(iday,:)],channels,[],N2Skip);
            end

            h=htmp;
            d = [d dtmp];

        end

    end
    if isempty(d),
        %disp(['no data read from the file of day ' num2str(iday)]);
        %disp(files(iday,:));
        %if isempty(iday),$disp('iday is empty');end
        %pause
    end


else

    d=[];
    h=[];

    disp(['File Not Found: ']);

end


fclose('all');


return

