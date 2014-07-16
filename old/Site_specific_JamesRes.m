function [D] = Site_specific_JamesRes(HEADER, D)

%This is a site specific file for James Reserve

%use time to remove bad data
time = D(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%By measurement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------
%Tsonic
%----------------------------------------------------------------------
%clean mat
%clean dl
%clean GOES
bad_Ts_goes=D(210,:)==0;
D(210,bad_Ts_goes)=NaN;