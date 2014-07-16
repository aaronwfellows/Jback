function [new,filt_arr]=multifilter(old,nsd,mean_or_med,Vmin,Vmax,disp_option,remove_hard_zero)
%function [new,filt_array]=multifilter(old,nsd,mean_or_med,Vmin,Vmax);
%Inputs
%%      old: a nx1 array that contains unfiltered data and possible NaNs
%%      nsd: number of std devns that define an outlier, enter 'x' if you
%%              don't want a nsd filter
%%      mean_or_med:  Enter 0 to use mean as measure of central tendency or 1 to use median
%%      Vmin: Force values to be no lower than this
%%      Vmax: Force values to be no higher than this
%%      disp_option (optional) Specify 1 if you want the results of
%%          filt_arr to be displayed (default is zero)
%%      remove_hard_zero option....specify 1 if you want the hard zeros
%%          removed
%Outputs
%%      new: the filtered array, returned in the same dimensions and length
%%              as old
%%      filt_array: a diagnostic array showing the number of each type of
%%      filter removed
%%
%%
%%      1. Length of original and newly filtered array
%%      2. Original number of NaNs
%%      3. New Number of NaNs (after filtering)
%%      4. Remaining non NaN values
%%      5. -99999
%%      6. -999
%%      7. +99999
%%      8. +9999
%%      9. number > Vmax
%%      10. Number < Vmin
%%      11. Number of nsd outliers
%%      12. Number of Hard Zeros removed
%       13. Percent of nonNan data removed by this filter
%%      
tmp=old;

if nargin<4,
    Vmin = -1e20;
    Vmax = 1e20;
    disp_option=0;
    remove_hard_zero=0;
end

if nargin<6,
    disp_option=0;
    remove_hard_zero=0;
end

if nargin<7,
    remove_hard_zero=0;
end

if ~isnumeric(nsd),
    nsd=-999;
end

arr_isnan=find(isnan(tmp));
num_isnan=length(arr_isnan);

arr_neg5nines=find(tmp==-99999);
num_neg5nines=length(arr_neg5nines);

arr_neg3nines=find(tmp==-999);
num_neg3nines=length(arr_neg3nines);

arr_pos5nines=find(tmp==-99999);
num_pos5nines=length(arr_pos5nines);

arr_pos4nines=find(tmp==-9999);
num_pos4nines=length(arr_pos4nines);

arr_too_high=find(tmp>=Vmax);
num_too_high=length(arr_too_high);

arr_too_low=find(tmp<=Vmin);
num_too_low=length(arr_too_low);

if remove_hard_zero==1,
    arr_zero=find(tmp==0);
    num_zero=length(arr_zero);
    tmp(arr_zero)=NaN;
else
    num_zero=NaN;
end

%replace all these with NaN and then run a nstdev filter over it

tmp(arr_neg5nines)=NaN;
tmp(arr_neg3nines)=NaN;
tmp(arr_pos5nines)=NaN;
tmp(arr_pos4nines)=NaN;
tmp(arr_too_high)=NaN;
tmp(arr_too_low)=NaN;


if nsd~=-999,
    if mean_or_med==0
        centtend=mnn(tmp);
    else
        centtend=nanmedian(tmp);
    end


    tmpdev=abs(tmp-centtend);
    critdev=nsd*nanstd(tmp);
    arr_ol=find(tmpdev>critdev);
    num_ol=length(arr_ol);

    tmp(arr_ol)=NaN;
else
    num_ol=NaN;
    if disp_option==1,
        disp('No Stdev outlier filtering');
        

    end


end
new=tmp;


%      1. Length of original and newly filtered array
%%      2. Original number of NaNs
%%      3. New Number of NaNs (after filtering)
%%      4. Remaining non NaN values
%%      5. -99999
%%      6. -999
%%      7. +99999
%%      8. +9999
%%      9. number > Vmax
%%      10. Number < Vmin
%%      11. Number of nsd outliers
%%      12. Percent of nonNan data removed by this filter
filt_arr(1:14)=NaN;

filt_arr(1)=length(new);
filt_arr(2)=num_isnan;
filt_arr(3)=length(find(isnan(new)));
filt_arr(4)=lnn(new);
filt_arr(5)=num_neg5nines;
filt_arr(6)=num_neg3nines;
filt_arr(7)=num_pos5nines;
filt_arr(8)=num_pos4nines;
filt_arr(9)=num_too_high;
filt_arr(10)=num_too_low;
if Vmin>-999,filt_arr(10)=num_too_low-num_neg3nines;end
filt_arr(11)=num_ol;
filt_arr(12)=num_zero;
filt_arr(13)=length(find(isnan(new)))-num_isnan;
filt_arr(14)=  100* (length(find(isnan(new)))-num_isnan) / (length(new)-num_isnan);

if disp_option==1
    FilterResults=strvcat('1. Length of original and newly filtered array',...
        '2. Original number of NaNs',...
        '3. New Number of NaNs (after filtering)',...
        '4. Remaining non NaN values',...
        '5. -99999',...
        '6. -999',...
        '7. +99999',...
        '8. +9999',...
        '9. number > Vmax',...
        '10. Number < Vmin',...
        '11. Number of nsd outliers',...
        '12. Number of hard zeros',...
        '13. Number of elements removed by this filter',...
        '14. Percent of nonNan data removed by this filter');
    disp('-------------------------------------------------');
    disp('------       MultiFilter Results             ----');
    disp('-------------------------------------------------');
    disp(['  # stdevs= ' num2str(nsd) ' , ' ' Vmin = ' num2str(Vmin) ' , ' ' Vmax = ' num2str(Vmax)]);
    disp('-------------------------------------------------');
    for ifx=1:14,


        disp([FilterResults(ifx,:) ' = ' num2str(filt_arr(ifx))]);

    end
    disp('-------------------------------------------------');
end

