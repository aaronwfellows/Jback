function [FLAG]=despike(x,nstds,xmin,xmax,signal)
%function [FLAG]=despike(x,nstds,xmin,xmax,signal)
%
% FLAG is a vector of 1's for values that do not
% qualify as spikes and 0's for values that do 
% qualify as spikes
%
%disp('Need to more carefully despike LI7500');

FLAG = ones(size(x));

idum = find(isnan(x));
FLAG(idum) = zeros(size( idum ) );

%aek if length(find(FLAG))>0
if ~isempty(find(FLAG, 1))    
   x(idum) = median(x(find(FLAG)))*ones(size(idum)); %#ok<*FNDSB> aek
   nspikestot=length(idum);
    if nspikestot > 0
        if nargin > 4
                fprintf('%-20s%8.0f',['# OF SPIKES [' signal ']: '],nspikestot);
        else
                fprintf('%-20s%8.0f','# OF SPIKES FOUND:',nspikestot);
        end
    end
else
   return
end

idum=find(x>xmax);    
FLAG( idum ) = zeros(size( idum ) );

%aek if length(find(FLAG))>0
if ~isempty(find(FLAG, 1))
   x(idum) = median(x(find(FLAG)))*ones(size(idum));
   nspikestot=nspikestot+length(idum);
   if nspikestot > 0
       fprintf('%8.0f','# of spikes found 2: ',nspikestot);
   end
else
   return
end

idum=find(x<xmin);    
FLAG( idum ) = zeros(size( idum ) );

%if length(find(FLAG))>0
if ~isempty(find(FLAG, 1))
    x(idum) = median(x(find(FLAG)))*ones(size(idum));
   nspikestot=nspikestot+length(idum);
   if nspikestot > 0; fprintf('%8.0f','# of spikes found 3: ',nspikestot); end;
else
   return
end


nspikes = 1;

while nspikes > 0
   
   idum=find( abs(x - median(x(find(FLAG)))) > nstds * std(x(find(FLAG))) );
   FLAG( idum ) = zeros(size( idum ) );    
   
   if ~isempty(find(FLAG, 1))
      x(idum) = median(x(find(FLAG)))*ones(size(idum));
   else
      fprintf('\n');
      return
   end
   
   nspikes=length(idum);
   
   if nspikes > 0
      nspikestot=nspikestot+nspikes;
      fprintf('%8.0f','# of spikes found 3: ',nspikestot);
   end
   
   %subplot(211);
   %plot(x(find(~FLAG)),'.'  )
   %hold on;
   %plot([1 length(find(~FLAG))],[median(x(find(FLAG)))+6*std(x(find(FLAG))) median(x(find(FLAG)))+6*std(x(find(FLAG)))],'r')
   %plot([1 length(find(~FLAG))],[median(x(find(FLAG)))-6*std(x(find(FLAG))) median(x(find(FLAG)))-6*std(x(find(FLAG)))],'r')
   %hold off
   
   %subplot(212);
   %plot(x(find(FLAG)),'.'  )
   %hold on;
   %plot([1 length(x(find(FLAG)))],[median(x(find(FLAG)))+6*std(x(find(FLAG))) median(x(find(FLAG)))+6*std(x(find(FLAG)))],'r')
   %plot([1 length(x(find(FLAG)))],[median(x(find(FLAG)))-6*std(x(find(FLAG))) median(x(find(FLAG)))-6*std(x(find(FLAG)))],'r')
   %hold off
   
   
   
end

fprintf('\n');


return






