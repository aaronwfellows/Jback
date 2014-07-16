function [header,data]=read5000(file,channels,Samples,N2Skip)

% July 21, 2002 -  modify to read partial file

if nargin < 4
    N2Skip=0;
end

if nargin < 3
    Samples = [];
end

if nargin < 2
    channels=[];
end

%CR=13;
LF=10;
COMMA=44;

fid=fopen(file,'r','ieee-le');
%disp(['fid = ' num2str(fid) '   file to open is ' file]);
d=fread(fid,10000,'uchar');

% find Line Feeds and Carriage Returns
%icr = find(d==CR);
ilf = find(d==LF);

% find the end of the header
if length(ilf)<5,
    %disp('describe ilf');
    whos ilf
end
EOH = ilf(5);

% read last line of header to get the file structure
HLine2 = d(ilf(1)+1:ilf(2)-1)';
HLine5 = d(ilf(4)+1:ilf(5)-1)';

begfields = [1 find(HLine5==COMMA)+1];
endfields = [find(HLine5==COMMA)-1 length(HLine5)-1];

begfields2 = [1 find(HLine2==COMMA)+1];
endfields2 = [find(HLine2==COMMA)-1 length(HLine2)-1];


Nfields = length(begfields);

for i=1:Nfields
    % don't read the quotes at beginning and end of each field
    FieldName{i} = char(HLine2(begfields2(i)+1:endfields2(i)-1));
    Field{i} = char(HLine5(begfields(i)+1: endfields(i)-1));
end


% Calculate the number of bytes in a record and get the
% corresponding matlab precision

for i=1:size(Field,2)
    if strcmp(char(Field(i)),'ULONG')
        NBytes(i) = 4; %#ok<*AGROW>
        MatlabPrec{i}='uint32';
    elseif  strcmp(char(Field(i)),'IEEE4L')
        NBytes(i) = 4;
        MatlabPrec{i}='float32';
    end
end

%% Start reading the channels

header = [];

% first position pointer at the end of the header
fseek(fid,EOH,'bof');
ftell(fid);

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

fprintf(1,'Reading %4.0f channels from [ %s ]',size(channels,1),file);
fprintf(1,'\n');

%%% if size(channels,1) >= 10
%%%    fprintf(1,'\n');
%%% end

NColumns=30+length(file);

for j=1:size(channels,1)
    %disp(['from line 135 of read5000 ' num2str(size(channels))]);
    ich=strmatch(channels(j,:),FieldName,'exact');

    Beg = EOH + BytesCumulative(ich) + BytesPerRecord*(Samples(1)-1);
    Beg=unique(Beg);
    fseek(fid,Beg,'bof');
    data(:,j)= fread(fid,N2Read,char(MatlabPrec(ich)),(N2Skip+1)*BytesPerRecord-NBytes(ich)) ; %#ok<AGROW>
    header = strvcat(header, char(FieldName{ich})); %#ok<VCAT>
    
    %%% fprintf(1,' .');
    NColumns=NColumns-2;
    if NColumns <= 0
        NColumns=30+length(file);
    end

end
%%% fprintf(1,'\n');
data=data';
fclose(fid);

return
