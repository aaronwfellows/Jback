function [dchan,index] = chset(d,ch,chstring,key)
%function [dchan,index] = chset(d,ch,chstring,key);
%
% 10/10;2001 - added optinal argument to return only the channel
% index and not the signal...
% key - 'index'  --> returns only the row number
% key = 'sig'    --> returns signal (default)

% given an MXN data array where M is the number of channels
% and N the number of samples, and ch an MX7 string array with
% the channel names corresponding to d, this function
% finds the location of channel 'chstring' in the ch array, and
% creates a vector (1XN) with the corresponding data.

if nargin<4
    key='sig';
end

casesensitive=0;

while length(chstring) < size(ch,2)
    chstring=[chstring ' ']; %#ok<*AGROW>
end

for idum=1:size(ch,1);
    if casesensitive==1,
        a(idum)=strcmp(ch(idum,:),chstring);
    else
        a(idum)=strcmpi(ch(idum,:),chstring);
    end
end

if findstr(key,'index')
    dchan=find(a);
else
    index=find(a);
    if isempty(index),
        disp(['Variable named ' chstring ' not found in chset....paused']);
        pause
    else
        dchan=d(index,:);
    end

end

