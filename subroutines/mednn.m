function [median_isnotanan]=mednn(X)

if ~isempty(find(~isnan(X), 1))
%if length(find(~isnan(X)))~=0,
   median_isnotanan=median(X(~isnan(X)));
else
   median_isnotanan=NaN;
end
return

