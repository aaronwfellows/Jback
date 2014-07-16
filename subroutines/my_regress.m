function [b,bint,r,rint,stats] = my_regress(y,x)

% b(1) is the intercept, b(2) is the slope
% stats = [r2 F prob s2]


if length(x)~=size(x,1),
    x=x';
end

if length(y)~=size(y,1),
    y=y';
end


%add a vector of ones adjacent to x

X = [ones(length(x),1) x];

[b,bint,r,rint,stats]=regress(y,X);
% 
% intercept = b(1)
% slope = b(2)

% 95% CI on intercept = bint(1, 1:2)
% 95% CI on slope = bint(2, 1:2)

% r = residuals

% rint = n x 2 array , a 95% CI for each residual

% stats = [r2 F prob s2]