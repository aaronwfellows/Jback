function [len_isnotanan]=lnn(X)

len_isnotanan=length(find(~isnan(X)));
