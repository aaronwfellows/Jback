function [h,d] = concat5000(Directory,channels,interval,N2Skip)

if nargin < 4
    N2Skip=0;
end

MaxSecOff = interval(3);

[FileBeg,SampleBeg]=FindSample5000(Directory,interval(1),MaxSecOff);
[FileEnd,SampleEnd]=FindSample5000(Directory,interval(2),MaxSecOff);
% disp(['interval is ' num2str(interval)])
if (~isempty(findstr('ts_data',Directory)) && (isempty(SampleBeg) || isempty(SampleEnd)))
    disp('Check FindSample5000 program, missing ts_data sampleBeg & sampleEnd');
    disp(Directory);
    %pause;
elseif (~isempty(findstr('flux',Directory))  && (isempty(SampleBeg) || isempty(SampleEnd)))
    disp('Check FindSample5000 program, missing flux sampleBeg & sampleEnd')
    disp(Directory);
    %pause;
elseif (~isempty(findstr('variance',Directory))  && (isempty(SampleBeg) || isempty(SampleEnd)))
    disp('Check FindSample5000 program, missing variance sampleBeg & sampleEnd')
    disp(Directory);
%%% What's the point of this part?
% else
%     disp(Directory);
%     disp('concat5000: Directory does not match file type in Findsample')
%     pause;
    %pause;
end

if (~isempty(FileBeg) && ~isempty(SampleBeg) && ~isempty(FileEnd) && ~isempty(SampleEnd))

    %  get a list of the data files in the directory
    %put in some logic to avoid unix measuring the width of the files array
    %incorrectly
    if ~isempty(findstr('ts_data',Directory)),
        fnamelength=27;
    elseif ~isempty(findstr('flux',Directory)),
        fnamelength=24;
    elseif ~isempty(findstr('variance',Directory)),
        fnamelength=28;
    else
        disp(['Directory does not match file type in Findsample' FileBeg Directory]);
        pause;
    end

    %now reduce the width of the files array to the length of the expected
    %filename
    files = dir(Directory);
    files = char(files.name);
    files=files(:,1:fnamelength);
    files_lower_names = lower(files);%this variable created for .tob comparison

    itablefile=strmatch('.tob',files_lower_names(:,end-3:end));

    files = files(itablefile,:);

    % sort the files by date
    filedates = str2num(files(:,[1:4 6 7 9 10 12:15])); %#ok<ST2NM>
    [WonkyStuff,index] = sort(filedates);

    files = files(index,:);

    ifilebeg=strmatch(FileBeg,files);
    ifileend=strmatch(FileEnd,files);

    d=[];
    h=[];

    if ifilebeg==ifileend
       [h,d] = read5000([Directory files(ifilebeg,:)],channels,[SampleBeg SampleEnd],N2Skip);
    else
        %disp('ifilebeg<>ifileend');
        for iday=ifilebeg:ifileend
            %disp('now reading the complete required data set, we are at l62 of concat5000 about to gointo read5000 ');
            if iday==ifilebeg
                [htmp,dtmp] = read5000([Directory files(iday,:)],channels,[SampleBeg NaN],N2Skip);
            elseif iday==ifileend
                [htmp,dtmp] = read5000([Directory files(iday,:)],channels,[1 SampleEnd],N2Skip);
            else
                [htmp,dtmp] = read5000([Directory files(iday,:)],channels,[],N2Skip);
            end

            h=htmp;
            d = [d dtmp]; %#ok<AGROW>
        end
    end
    if isempty(d)
        disp(['no data read from the file of day ' num2str(iday) ' in concat5000']);
        pause;
    end
else
    d=[];
    h=[];
    %disp(['File Not Found in concat5000: ' FileBeg FileEnd Directory interval]);
    %pause;
end


fclose('all');


return

