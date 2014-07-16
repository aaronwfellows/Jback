function [] = Rn_patch_Sage()
%Make a new Rn for Sage

%The sage has a lower quality net radiation (Rn) sensor.  We are going to
%use the Grassland Rn to construct a time series of Rn over the sage.  The
%biggest difference between the sites is the albedo.  We need to adjust the 
%outgoing shortwave radiation.

%------------------------------------------------------------------------%
%Grab the data
load('C:\towerData\WebsiteData\Sage_v3_2.mat');

CSS=DATA;
H_CSS=HEADER;

load('C:\towerData\WebsiteData\Grass_v3_2.mat');

GRASS=DATA;
H_GRASS=HEADER;

%------------------------------------------------------------------------%
%Radiation data
len=size(CSS, 1);

%[CSS time stamp; Grass time stamp; Rn Sage; Rn Grass; SW out Sage; SW out_Grass; SW in Sage; SW in_Grass; Rn fill];
RadiationHeader = {'CSS_time'; 'Grass_time'; 'Rn_CSS'; 'Rn Grass'; 'SW_out_CSS'; 'SW out_Grass'; 'SW_in_CSS'; 'SW_in_Grass'; 'Rn_fixed_used_in_web'};
Radiation = ones(len,9) *NaN;
%------------------------------------------------------------------------%
%The datasets are not lined up - Merge
%------------------------------------------------------------------------%
%We want to get to a new CSS Rn.  Use the CSS time stamp to line up the
%data
%time
CSS_time = CSS(:,1);
GRASS_time = GRASS(:,1);

%convert tiem to integer
CSS_int = round(CSS_time*48);
GRASS_int = round(GRASS_time*48);

[C,i_css,i_grass] = intersect(CSS_int,GRASS_int);%C is the value in common, i_css is the column in css that 

%common radiation
Radiation(:,1) = CSS(:,1);%time
Radiation(i_css,2) = GRASS(i_grass,1);%time
Radiation(:,3) = CSS(:,15);%Rn
Radiation(i_css,4) = GRASS(i_grass,15);%Rn
Radiation(:,5) = CSS(:,14);%out
Radiation(i_css,6) = GRASS(i_grass,14);%out
Radiation(:,7) = CSS(:,13);%in
Radiation(i_css,8) = GRASS(i_grass,13);%in
%------------------------------------------------------------------------%
%albedo for fill
%albedo = outgoing SW/incoming SW

%Grass
good = Radiation(:,6) > 50 & Radiation(:,8) > 300;
albedo_g = Radiation(good,6)./Radiation(good,8);
albedo_grass = nanmedian(albedo_g)

%Use this to make a filled SW in and SW out 
%This does not work
%miss_6=isnan(Radiation(:,6));

%fill_6 = Radiation(:,8) * albedo_grass;
%Radiation(miss_6,6) = fill_6(miss_6,1);


%Sage - this works better
%outgoing radiation from grass is probably OK for the periods that are
%missing - the largest missing period follows the burn when it is likely
%that herbaceous annual dominanted the veg - while this is not perfect it
%is likely better than other options
miss = isnan(Radiation(:,5))==1 | isnan(Radiation(:,6))==1;
P = polyfit(Radiation(~miss,5),Radiation(~miss,6),1); %P = polyfit(X,Y,N)
fill_5 = (Radiation(:,6)- P(1,2))/ P(1,1); 


miss_5=isnan(Radiation(:,5));
Radiation(miss_5,5) = fill_5(miss_5,1);
%------------------------------------------------------------------------%
%fix up the Rn
Rn_fix = Radiation(:,4) + Radiation(:,6) - Radiation(:,5);

Radiation(:,9) = Rn_fix(:,1);
%------------------------------------------------------------------------%
%compare Rn_fix with Rn at the CSS

figure(1)
plot(CSS(:,1), CSS(:,15), 'k.')
hold on
plot(GRASS(:,1), GRASS(:,15), 'bo')
hold on
plot(Radiation(:,1), Radiation(:,9), 'rx')
%------------------------------------------------------------------------%
%------------------------------------------------------------------------%
%save radiation
site_name = 'Sage'
out_root='C:\towerData\Rn_fix\';
filename = [out_root site_name '_' 'Radiation' '_' 'v3' '_' '2']

DDD=date;

%save as a mat file
save(filename, 'Radiation', 'RadiationHeader', 'DDD');

%------------------------------------------------------------------------%
%------------------------------------------------------------------------%
%Patch up CSS DATA file
CSS(:,15) = Radiation(:,9);
DATA = CSS;
HEADER = H_CSS;
%------------------------------------------------------------------------%
%save 'C:\towerData\WebsiteData\Sage_v3_2.mat'
site_name = 'Sage'
out_root='C:\towerData\WebsiteData\';
filename = [out_root site_name '_' 'v3' '_' '2']

%save as a mat file
save(filename, 'DATA', 'HEADER');
%------------------------------------------------------------------------%
%------------------------------------------------------------------------%