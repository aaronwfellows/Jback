function [xnorm,ynorm,slope,offset,r2]=normalize_filt(x_orig,y_orig,plotoption,covarmaxoption,regfilteroption)
% function  [xnorm,ynorm,slope,offset,r2]=normalize(x_orig,y_orig,plotoption,covarmaxoption,regfilteroption);
%
%
% adjust y so that it is normalized to x and return the linear coefficients
% to transform the original y to the new y
%
%
% INPUTS:
%           x_orig: a n x 1 array of data that will NOT be adjusted
%           y_orig: a n x 1 array, correlated to x_orig, but different by
%                       an offset and a gain which this function will
%                       estimate
%           plotoption:
%                set to 1 if you want to plot the data and inspect the
%                adjustment, set to 2 if you also want to see the effects
%                of reg filter
%
%           covarmaxoption: set this to one if you suspect that there is a
%               delay between the arrays
%
%           regfilteroption: set this to 1 if you want to filter outliers
%               from the regression curve
%
% OUTPUTS:
%           xnorm: same as x_orig
%           ynorm: normalized data (same dimensions as input
%           slope,offset r2: linear relationship between x and y
%
%
%
%
%

%  Functions called: covarmax2.m, regfilter.m, myreg.m
%
%
%close all;
%[x,y]=eatnan(x_orig,y_orig);
x=x_orig;
y=y_orig;

goodbix = ~isnan(x) & ~isnan(y) & x~=-999 & y~=-999 & x~=-972.3727 & y~=-972.3727 & x~=0 & y~=0;

x=x(goodbix);
y=y(goodbix);

if covarmaxoption==1,
    %[cmax,shift,x,y]=covarmax2(x,y,10);
    [~,~,x,y]=covarmax2(x,y,10);
end

if regfilteroption==1,

    if plotoption==2,
        plotflag=1;
    else
        plotflag=0;
    end

    [y]=regfilter(x,y,2,plotflag);
end



[m,b,Sy,Sm,Sb,r2]=myreg(x,y);

yfit=m*x+b;

if plotoption==1 || plotoption==2 
    figure(3001);
    plot(x,y,'k.');
    hold on;
    plot(x,yfit,'r');

    ctext(90,10,['Slope = ' num2str(m) ', Intercept = ' num2str(b) ', r2 = ' num2str(r2)]);
    ctext(80,10,'Hit return to plot normalized data');
    hold off
    pause
    clf;
    close(3001);
end

yadj = (y - b) /m;

[m_2,b_2,Sy_2,Sm_2,Sb_2,r2_2]=myreg(x,yadj);
yfit_2=m_2*x+b_2;


if plotoption==1 || plotoption==2 ,
    figure(3002);
    H=plot(x,yadj,'b.');
    hold on;
    plot(x,yfit_2,'m','linewidth',3);




    Lims= [ get(gca,'XLim') ; get(gca,'YLim') ];
    newminlim = min(Lims(:,1));
    newmaxlim = max(Lims(:,2));
    newlim = [newminlim newmaxlim];
    limrange = newmaxlim-newminlim;
    nsteps = 50;
    pnts = limrange / nsteps;

    vec = newminlim:pnts:newmaxlim;

    plot(vec,vec,'k:','linewidth',2);


    set(gca,'YLim',newlim,'XLim',newlim);

    ctext(95,10,'Normalized Data');
    ctext(90,10,['Slope = ' num2str(m_2) ', Intercept = ' num2str(b_2) ', r2 = ' num2str(r2_2)]);
    ctext(85,10,['Old Slope = ' num2str(m) ', Old Intercept = ' num2str(b) ]);
    ctext(80,10,'Hit return to finish');
    hold off
    pause
end

slope = m;
offset = b;
%r2 = r2;
xnorm = x_orig;
ynorm = (y_orig - offset) / slope;


chksection=0;

if chksection==1,
    x1=x_orig;
    y1=y_orig;

    len = length(x1);
    nints = 10;
    step = fix(len/nints);

    for i=1:nints,
        ix =  (i-1) * step + 1 : i * step ;
        xc=x1(ix);
        yc=y1(ix);

        figure(999);
        plot(xc,yc,'r.');hold on
        plot(xc,xc,'g');hold off
        
        minx = min([      min(xc(xc>-40))        min(yc(yc>-40))         ]);
        maxx = max([     max(xc(xc<-40))         max(yc(yc<40))         ]);
        
        axis([minx maxx minx maxx]);
        
        
        
        ctext(90,10,[num2str(min(ix)) ' to ' num2str(max(ix))    ]  );
        pause
    end
end






return
