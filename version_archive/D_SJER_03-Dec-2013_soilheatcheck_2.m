function [Q_Rn_polyfitRegression] = soilheatcheck_2(Cleaned_D, Ts, delTs, Vwater, Vs, Cvs, Qs, delQ, site_name, out_root)

%------------------------------------------------------------------------%
%raw time series plots
%------------------------------------------------------------------------%
time = Cleaned_D(:,1);

figure(1)
plot(time(:,1), Ts(:,1), 'r.')
hold on
plot(time(:,1), Ts(:,2), 'y.')
hold on
plot(time(:,1), Ts(:,3), 'm.')
hold on
plot(time(:,1), Ts(:,4), 'g.')
hold on
plot(time(:,1), Ts(:,5), 'b.')
ylabel('Ts (C)')
xlabel('time')
legend('30', '34', '50', '100', '200')
title(site_name)

f='Ts_timeseries';
saveas(figure(1), [out_root f]);

figure(2)
plot(time(:,1), delTs(:,1), 'r.')
hold on
plot(time(:,1), delTs(:,2), 'y.')
hold on
plot(time(:,1), delTs(:,3), 'm.')
hold on
plot(time(:,1), delTs(:,4), 'g.')
hold on
plot(time(:,1), delTs(:,5), 'b.')
ylabel('delTs (C)')
xlabel('time')
legend('30', '34', '50', '100', '200')
title(site_name)

f='delTs_timeseries';
saveas(figure(2), [out_root f]);

figure(3)
plot(time(:,1), Vs(:,1), 'r.')
hold on
plot(time(:,1), Vs(:,2), 'b.')
hold on
plot(time(:,1), Vs(:,3), 'm.')
hold on
plot(time(:,1), Vs(:,4), 'g.')
hold on
plot(time(:,1), Vwater(:,1), 'k.')
ylabel('V of water (fraction)')
xlabel('time')
legend('1', '2', '3', '4', 'nanmean')
title(site_name)

f='Vwater_timeseries';
saveas(figure(3), [out_root f]);

figure(4)
plot(time(:,1), Cvs(:,1), 'ro')
hold on
plot(time(:,1), Cvs(:,2), 'yx')
hold on
plot(time(:,1), Cvs(:,3), 'm.')
hold on
plot(time(:,1), Cvs(:,4), 'go')
hold on
plot(time(:,1), Cvs(:,5), 'bx')
ylabel('Cvs (J/C/m^3)')
xlabel('time')
legend('30', '34', '50', '100', '200')
title(site_name)

f='Cvs_timeseries';
saveas(figure(4), [out_root f]);

figure(5)
plot(time(:,1), Qs(:,1), 'r.')
hold on
plot(time(:,1), Qs(:,2), 'y.')
hold on
plot(time(:,1), Qs(:,3), 'm.')
hold on
plot(time(:,1), Qs(:,4), 'g.')
hold on
plot(time(:,1), Qs(:,5), 'b.')
ylabel('Qs (W/m^2)')
xlabel('time')
legend('30', '34', '50', '100', '200')
title(site_name)

f='Qs_timeseries';
saveas(figure(5), [out_root f]);

figure(6)
plot(time(:,1), delQ(:,1), 'g.')
hold on
plot(time(:,1), Cleaned_D(:,15), 'k.')
hold on
plot(time(:,1), Cleaned_D(:,10), 'b.')
hold on
plot(time(:,1), Cleaned_D(:,9), 'r.')
ylabel('(W/m^2)')
xlabel('time')
legend('Q soil storage', 'Rn', 'LE', 'SH')
title(site_name)

f='EnergyFluxes';
saveas(figure(6), [out_root f]);


figure(7)
plot(time(:,1), delQ(:,1)./Cleaned_D(:,15), 'k.')
ylabel('Q soil / Rn ')
xlabel('time')
legend('Q soil storage/Rn')
title(site_name)

f='EnergyFluxes';
saveas(figure(7), [out_root f]);


%find how much of Rn goes to soil heat storage
usable_obs = isfinite(Cleaned_D(:,15)) ==1 & isfinite(delQ(:,1)) ==1;
P=polyfit(Cleaned_D(usable_obs,15), delQ(usable_obs,1), 1)
Q_Rn_polyfitRegression = P;

figure(8)
plot(Cleaned_D(:,15), delQ(:,1), 'k.')
hold on
plot(Cleaned_D(:,15), P(1,1) * Cleaned_D(:,15) + P(1,2), 'r-')
ylabel('Q soil')
xlabel('Rn')
title(site_name)

f='EnergyFluxes';
saveas(figure(8), [out_root f]);

pause
%------------------------------------------------------------------------%
%The daily cycle
%------------------------------------------------------------------------%
%GMT to California time
Caltime = time - (8/24); %[days]
day = floor(Caltime);
f_day = (Caltime - day); %[frac day]

%30 min time step
inc=1/48;

%turbulent periods
turb=Cleaned_D(:,4) > 0.3 & isfinite(Cleaned_D(:,9)) & isfinite(Cleaned_D(:,10)) & isfinite(Cleaned_D(:,15)) & isfinite(delQ(:,1));

%Daily heat flux patterns
Daily_Heat =ones(48,5)*NaN;

for i=1:48
    t=(i*inc)-inc;
    good = f_day > (t - 0.001) & f_day < (t + 0.001) & turb(:,1) == 1;
    d_G=mean(delQ(good,1));
    d_Rn=mean(Cleaned_D(good,15));
    d_H=mean(Cleaned_D(good,9));
    d_LE=mean(Cleaned_D(good,10));
    Daily_Heat(i,1)=t;
    Daily_Heat(i,2)=d_G;
    Daily_Heat(i,3)=d_Rn;
    Daily_Heat(i,4)=d_H;
    Daily_Heat(i,5)=d_LE;
end

figure(9)
plot(Daily_Heat(:,1), Daily_Heat(:,2), 'g.')
hold on
plot(Daily_Heat(:,1), Daily_Heat(:,3), 'k.')
hold on
plot(Daily_Heat(:,1), Daily_Heat(:,4), 'r.')
hold on
plot(Daily_Heat(:,1), Daily_Heat(:,5), 'b.')
ylabel('(W/m^2)')
xlabel('time')
legend('Q soil storage', 'Rn', 'SH', 'LE')
title(site_name)
 
f='DailyCycle_EnergyFluxes';
saveas(figure(9), [out_root f]);

pause

close(figure(1))
close(figure(2))
close(figure(3))
close(figure(4))
close(figure(5))
close(figure(6))
close(figure(7))

copyfile('soilheatcheck_2.m',['./version_archive/' site_name '_' date '_soilheatcheck_2.m'],'f');



