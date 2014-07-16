function [cmax,N] = covarmax(a,b,ShiftLim,P) 

% Find sample shift for  maximum covariance between A and B.
% Two Cases are considered, depending on whether ShiftLim has
% one or 2 elements:
%
%     1) Case1: Shiftlim has 2 elements: [minshift maxshift].  
%        It is known that signal B lags A in time, so B is shifted 
%        for sample lags from minshift to maxshift. The number N is
%        the sample shift of B w.r.t. A
%
%     2) Case 2: ShiftLim has 1 elment - the maximum expected shift
%        for either B or A.  Here, A is shifted relative to B, then
%        B is shifted relative to A.  N is a 2 element vector, the
%        first element is the shift of A relative to B, and the 
%        second is the shift of B relative to A.  One of the elements
%        will be zero, unless there was a tie for maximum covariance,
%        which I'm not sure how to interpret
%

   %%%%%%%%%%%%%%%%%%%%%%%
   %%% CASE 1: B lags A
   %%%%%%%%%%%%%%%%%%%%%%%


   if length(ShiftLim)==2

     [c,cmax,N] = docov(a,b,(ShiftLim(1):1:ShiftLim(2)));      
     ShiftAx = (ShiftLim(1):1:ShiftLim(2));
     Comment='B Lags A';

   elseif length(ShiftLim)==1

   %%%%%%%%%%%%%%%%%%%%%%%
   %%% CASE 1: B lags A
   %%%%%%%%%%%%%%%%%%%%%%%

     [ca,camax,Na] = docov(a,b,(0:1:ShiftLim));      
     [cb,cbmax,Nb] = docov(b,a,(0:1:ShiftLim));      

     if abs(camax)>abs(cbmax)

        N=(Na);
        c=ca;
        ShiftAx = (0:1:ShiftLim);
        Comment='B Lags A';
        cmax=camax;
 
     elseif abs(cbmax)>abs(camax)

       N=(-Nb);
       c=cb;
       ShiftAx = (0:1:ShiftLim);
       Comment='A Lags B';
       cmax=cbmax;

     elseif camax==cbmax

      N=(Na);
      c=[fliplr(cb) ca];
      ShiftAx = [-(length(cb)-1):0 0:length(cb)-1];
      Comment='A and B tied';
      cmax=camax;

     end

   end


      if nargin > 3


     %%FIGURE subplot

     eval([ 'figure(' int2str(P.figure) ');']);
     eval([ 'subplot' P.sub ';']);

     plot(ShiftAx,c,'.');
     set(gca,'fontsize',8);
     title( [ P.title '  ' Comment ':   '  int2str(N)]);
     %xlabel('ShiftAxis')
     ylabel('Cov')

     end
 


return


function [covab,cmax,N] = docov(sig1,sig2,Shift)

   len=length(sig2);

   for ilag=1:length(Shift)
    
     %  [1 len-Shift(ilag); Shift(ilag)+1 len];

       covi=cov(sig1(1:len-Shift(ilag)),sig2(Shift(ilag)+1:len));
       covab(ilag)=covi(1,2); %#ok<AGROW>
    
   end

       ishift = find(abs(covab)==max(abs(covab)), 1 );
       N=Shift(ishift);
       cmax = covab(ishift);
return
