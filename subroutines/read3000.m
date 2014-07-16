function [header,data]=read3000(file,channels,Samples,N2Skip)
%function [header,data]=read3000(file,channels,Samples,N2Skip)
%
% Oct 2006 - modifying to handle more header lines
%
% April 2006 - modifying to read ascii too
%
% reads CR3000 and CR5000 data, outputs a header and the data array
%
% written by: Scott Miller (2002-2007) 
%             Atmospheric Sciences Research Center
%             SUNY Albany
%
% inputs:
% file - full file path
% channels - optional input string array that specifies the
%            channels to be read (to extract data); for example
%            channels={'channel1';'channel2'}
% Samples - optional input to read a specific subset of the
%           data.  Input as [firstsample lastsample].
%           If Samples is empty, reads all of the data.
%           If lastsample is NaN, reads from firstsample to end
% N2Skip - optional input to skip N2Skip samples.  Does not
%          do averaging, simply skips samples
%
% Note: a useful use of this script it to call the function with
% no output variables specified - this will read the header and
% scroll it to the Matlab terminal window, and is useful to query
% what is in a file; i.e.,  >>read3000(file)
%
% April 2006: modified to read ascii files

if nargin < 4
  N2Skip=0;
end

if nargin < 3
  Samples = [];
end

if nargin < 2
  channels=[];
end

CR=13;
LF=10;
COMMA=44;

% added deblank Aug 9, 2006
fid=fopen(deblank(file),'r','ieee-le');
d=fread(fid,10000,'uchar');

% find Line Feeds and Carriage Returns

icr = find(d==CR);
ilf = find(d==LF);

%d(icr)

%size(icr)
%size(ilf)

% disp(['line 2: ' int2str(d([icr(2) icr(2)+1 icr(2)+2])')]);
% disp(['line 3: ' int2str(d([icr(3) icr(3)+1 icr(3)+2])')]);
% disp(['line 4: ' int2str(d([icr(4) icr(4)+1 icr(4)+2])')]);
% disp(['line 5: ' int2str(d([icr(5) icr(5)+1 icr(5)+2])')]);
% disp(['line 6: ' int2str(d([icr(6) icr(6)+1 icr(6)+2])')]);
% disp(['line 7: ' int2str(d([icr(7) icr(7)+1 icr(7)+2])')]);

switch d(4)
 case abs('A')

  EOH=ilf(4);
  
  % read last line of header to get the file structure
  HLine2 = d(ilf(1)+1:ilf(2)-1)';
  HLine3 = d(ilf(2)+1:ilf(3)-1)';
  HLine4 = d(ilf(3)+1:ilf(4))';

% $$$   disp(char(HLine2));
% $$$   disp(char(HLine3));
% $$$   disp(char(HLine4));

  begfields2 = [1 find(HLine2==COMMA)+1];
  endfields2 = [find(HLine2==COMMA)-1 length(HLine2)-1];

  begfields3 = [1 find(HLine3==COMMA)+1];
  endfields3 = [find(HLine3==COMMA)-1 length(HLine3)-1];

  begfields4 = [1 find(HLine4==COMMA)+1];
  endfields4 = [find(HLine4==COMMA)-1 length(HLine4)-1];

  Nfields = length(begfields2);

  for i=1:Nfields
    % don't read the quotes at beginning and end of each field
    FieldName{i} = char(HLine2(begfields2(i)+1:endfields2(i)-1));
    FieldUnits{i}= char(HLine3(begfields3(i)+1: endfields3(i)-1));
  end
  
 case abs('B')
  % nov 10, 2006: redoing this because it chokes when the first
  % character of the binary data after the last line of the header
  % turns out to be a double quote in ascii - this happens some
  % times!
  % 
  % find the end of the binary header - find the last occurrence of
  % the ascii combination [34 13 10] which is 
  % ['double quote' 'carriage return' 'line feed' ]  
  % this line is the last line of the header

  %  strfind(char(d'),char([34 13 10]))
  nrowsa=length(strfind(char(d'),char([34 13 10])));
  
  % double check the end of the binary header by looking for the ascii combination 
  % [34 13 10 34] which is 
  % ['double quote' 'carriage return' 'line feed' 'double quote' ]  
  % this line is the second last line of the header

  %  strfind(char(d'),char([34 13 10 34]))

  %disp('11-24-06: using the search for the second-last header row, instead of last row');
  nrowsb=length(strfind(char(d'),char([34 13 10 34])))+1;
  nrows=min([nrowsa nrowsb]);

  EOH=ilf(nrows);

% $$$   % find the end of the binary header - find the last occurrence of
% $$$   % the ascii combination [13 10 34] which is 
% $$$   % ['carriage return' 'line feed' 'double quote']  
% $$$   % The line following this line is the last line of the header
% $$$ 
% $$$   nrows=length(strfind(char(d'),char([13 10 34])))+1
% $$$ 
% $$$   EOH=ilf(nrows)
  
  % read last line of header to get the file structure
  HLineName = d(ilf(nrows-4)+1:ilf(nrows-3)-1)';
  HLineUnits = d(ilf(nrows-3)+1:ilf(nrows-2)-1)';
  HLinePrecision = d(ilf(nrows-1)+1:ilf(nrows))';

% $$$   disp(char(HLineName));
% $$$   disp(char(HLineUnits));
% $$$   disp(char(HLinePrecision));

  begfieldsName = [1 find(HLineName==COMMA)+1];
  endfieldsName = [find(HLineName==COMMA)-1 length(HLineName)-1];

  begfieldsUnits = [1 find(HLineUnits==COMMA)+1];
  endfieldsUnits = [find(HLineUnits==COMMA)-1 length(HLineUnits)-1];

  begfieldsPrecision = [1 find(HLinePrecision==COMMA)+1];
  endfieldsPrecision = [find(HLinePrecision==COMMA)-1 length(HLinePrecision)-2];

  Nfields = length(begfieldsName);

  for i=1:Nfields
    % don't read the quotes at beginning and end of each field
    FieldName{i}=char(HLineName(begfieldsName(i)+1:endfieldsName(i)-1));
    FieldUnits{i}=char(HLineUnits(begfieldsUnits(i)+1: endfieldsUnits(i)-1));
    FieldPrecision{i}=char(HLinePrecision(begfieldsPrecision(i)+1: endfieldsPrecision(i)-1));
  end
end

%EOHA = ilf(4);
%EOHB = ilf(5);

if nargout<1
  for i=1:Nfields
    switch d(4)
     case abs('A')
      disp([char(FieldName{i}) ', (' char(FieldUnits{i}) ')']);
     case abs('B')
      disp([char(FieldName{i}) ', (' char(FieldUnits{i}) '), ('  char(FieldPrecision{i}) ')']);
    end
  end
  return
end

% determine if file is ASCII or binary
switch d(4)
 case abs('A')
  header=[];
  for ifield=1:Nfields
    header=strvcat(header,char(FieldName{ifield}));
  end

  disp(['reading ASCII file ' file ' ...']);
  fseek(fid,0,'bof');
  ftell(fid);

  d=fread(fid,inf,'uchar');
  d=d(EOH+1:end); % remove header

  % find Line Feeds and Carriage Returns
  icr = find(d==CR);
  ilf = [0; find(d==LF)];

  nlines=length(ilf);

  data=cell(nlines,Nfields);

  % after the header, read each line of data
  for iline=1:nlines-1
    TextLine = d(ilf(iline)+1:ilf(iline+1)-1)';
    
    %     disp([int2str(iline) ': ' char(TextLine)]);
    begfields = [1 find(TextLine==COMMA)+1];
    endfields = [find(TextLine==COMMA)-1 length(TextLine)-1];

    NfieldsInLine = length(begfields);

    if NfieldsInLine==Nfields
      for ifield=1:Nfields
        % check for quotes - quotes read as string
        if isempty(findstr(char(TextLine(begfields(ifield): endfields(ifield))),'"'))
          %          disp('no quotes')
          data{iline,ifield} = char(TextLine(begfields(ifield): endfields(ifield)));
        else
          % don't read the quotes at beginning and end of each field
          data{iline,ifield} = char(TextLine(begfields(ifield)+1: endfields(ifield)-1));
        end
        %         disp(data(iline,ifield));
      end
    else
      %       disp(['Line ' int2str(iline) ' has ' int2str(NfieldsInLine) ' fields: needs to have ' int2str(Nfields)]);
    end
    %      pause

  end
  data=data';

 case abs('B')
  % Calculate the number of bytes in a record and get the
  % corresponding matlab precision
%  FieldPrecision
  for i=1:size(FieldPrecision,2)
    if strcmp(char(FieldPrecision(i)),'ULONG')
      NBytes(i) = 4;
      MatlabPrec{i}='uint32';
    elseif  strcmp(char(FieldPrecision(i)),'IEEE4L')
      NBytes(i) = 4;
      MatlabPrec{i}='float32';
    elseif  strcmp(char(FieldPrecision(i)),'IEEE4')
      NBytes(i) = 4;
      MatlabPrec{i}='float32';
    end
  end

  %%% Start reading the channels
  header = [];

  % first position pointer at the end of the header
  fseek(fid,EOH,'bof');ftell(fid);

  BytesPerRecord=sum(NBytes)*ones(size(NBytes)) ;
  BytesCumulative = [0 cumsum(NBytes(1:length(NBytes)-1))];

  % check if a subset of the data is to be read, or the entire file
  % should be read.  Or, if Samples is nonzero, but the second
  % element of Samples is NaN, then read to the end of the file, by
  % specifiying N2Read as inf
  if isempty(Samples)
    N2Read=inf;
    Samples=1;
  elseif isnan(Samples(2))
    N2Read = inf;
  else
    N2Read = (Samples(2)-Samples(1)+1)/(N2Skip+1);
  end

  if isempty(channels)
    channels = char(FieldName);
  end

  fprintf(1,'Reading %4.0f channels from [ %s ]:',size(channels,1),file);

  if size(channels,1) >= 10
    fprintf(1,'\n');
  end

  NColumns=30+length(file);

  for j=1:size(channels,1)

    ich=strmatch(deblank(channels(j,:)),FieldName,'exact');
    
    Beg = EOH + BytesCumulative(ich) + BytesPerRecord(1)*(Samples(1)-1);

    fseek(fid,Beg,'bof');

    %disp([int2str(j) ':  ' char(MatlabPrec(ich))]);

    data(:,j)= fread(fid,N2Read,char(MatlabPrec(ich)),(N2Skip+1)*BytesPerRecord-NBytes(ich)) ;

    %data(:,j)

    header=strvcat(header,char(FieldName{ich}));

    fprintf(1,' .');

    NColumns=NColumns-2;
    if NColumns <= 0
      fprintf(1,'\n');
      NColumns=30+length(file);
    end
  end

  fprintf(1,'\n');
  data=data';

end

fclose(fid);

return

