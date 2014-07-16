function [] = Rn_patch_PFBurn()
%Make a new Rn for PF_Burn

%Fix up the PF burn Rn.
%------------------------------------------------------------------------%
%Grab the data
load('..\WebsiteData\PinyonBurn_v3_2.mat');

PinyonBurn=DATA;
H_PinyonBurn=HEADER;

load('..\WebsiteData\Pinyon_v3_2.mat');

Pinyon=DATA;
H_Pinyon=HEADER;

%------------------------------------------------------------------------%
%Radiation data
len=size(PinyonBurn, 1);

RadiationHeader = {'PinyonBurn_time'; 'Pinyon_time'; 'Rn_PinyonBurn'; 'Rn Pinyon'; 'SW_out_PinyonBurn'; 'SW out_Pinyon'; 'SW_in_PinyonBurn'; 'SW_in_Pinyon'; 'Rn_fixed_used_in_web'};
Radiation = ones(len,9) *NaN;
%------------------------------------------------------------------------%
%The datasets are not lined up - Merge
%------------------------------------------------------------------------%
%We want to get to a new CSS Rn.  Use the CSS time stamp to line up the
%data
%time
PinyonBurn_time = PinyonBurn(:,1);
Pinyon_time = Pinyon(:,1);

%convert tiem to integer
PinyonBurn_int = round(PinyonBurn_time*48);
Pinyon_int = round(Pinyon_time*48);

[C,i_PinyonBurn,i_Pinyon] = intersect(PinyonBurn_int,Pinyon_int);%C is the value in common, i_css is the column in css that 

%common radiation
Radiation(:,1) = PinyonBurn(:,1);
Radiation(i_PinyonBurn,2) = Pinyon(i_Pinyon,1);
Radiation(:,3) = PinyonBurn(:,15);
Radiation(i_PinyonBurn,4) = Pinyon(i_Pinyon,15);
Radiation(:,5) = PinyonBurn(:,14);%out
Radiation(i_PinyonBurn,6) = Pinyon(i_Pinyon,14);%out
Radiation(:,7) = PinyonBurn(:,13);%in
Radiation(i_PinyonBurn,8) = Pinyon(i_Pinyon,13);%in
%------------------------------------------------------------------------%
%albedo for fill
%albedo = outgoing SW/incoming SW

%Pinyon
good = Radiation(:,6) > 50 & Radiation(:,8) > 300;
albedo_p = Radiation(good,6)./Radiation(good,8);
albedo_Pinyon = nanmedian(albedo_p);

%Use this to make a filled SW in and SW out 
miss_6=isnan(Radiation(:,6));

fill_6 = Radiation(:,8) * albedo_Pinyon;
Radiation(miss_6,6) = fill_6(miss_6,:);


%PinyonBurn
good = Radiation(:,5) > 50 & Radiation(:,7) > 300;
albedo_pb = Radiation(good,5)./Radiation(good,7); %PB out/PBin
albedo_PinyonBurn = nanmedian(albedo_pb); %mean
%Use this to fill SW out at Burn
fill_5 = Radiation(:,7) * albedo_PinyonBurn;

miss_5=isnan(Radiation(:,5));
Radiation(miss_5,5) = fill_5(miss_5,:);
%------------------------------------------------------------------------%
%fix up the Rn
Rn_fix = Radiation(:,4) + Radiation(:,6) - Radiation(:,5);

Radiation(:,9) = Rn_fix(:,1);
%------------------------------------------------------------------------%
%compare Rn_fix 

figure(1)
plot(PinyonBurn(:,1), PinyonBurn(:,15), 'ro')
hold on
plot(Pinyon(:,1), Pinyon(:,15), 'b.')
hold on
plot(Radiation(:,1), Radiation(:,9), 'kx')

%------------------------------------------------------------------------%
%save radiation
site_name = 'PinyonBurn'
out_root='..\Rn_fix\';
filename = [out_root site_name '_' 'Radiation' '_' 'v3' '_' '2']

DDD=date;

%save as a mat file
save(filename, 'Radiation', 'RadiationHeader', 'DDD');

%------------------------------------------------------------------------%
%Patch up PinyonBurn DATA file
PinyonBurn(:,15) = Radiation(:,9);
DATA = PinyonBurn;
HEADER = H_PinyonBurn;
%------------------------------------------------------------------------%
%save 'C:\towerData\WebsiteData\PinyonBurn_v3_2.mat'
site_name = 'PinyonBurn'
out_root='C:\towerData\WebsiteData\';
filename = [out_root site_name '_' 'v3' '_' '2']

%save as a mat file
save(filename, 'DATA', 'HEADER');
%------------------------------------------------------------------------%
%------------------------------------------------------------------------%