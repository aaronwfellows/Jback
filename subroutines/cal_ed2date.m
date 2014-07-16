function [cdate]=cal_ed2date(ExpDay)
% updated 9/9/10 aek
global towerYearStart iSite

cdate = towerYearStart(iSite) + ExpDay;
cdate=datestr(cdate,0);
return;
