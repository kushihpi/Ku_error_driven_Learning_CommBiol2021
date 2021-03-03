function  [p, qhat, sigsq, qhatold, sigsqold]=recfilter(I, sige, qguess,rhoone, delk, nuone, muone)

% Written by Anne C. Smith 2004
% publication: Smith, A. C. et al. Dynamic analysis of learning in behavioral experiments. The Journal of neuroscience : the official journal of the Society for Neuroscience 24, 447-461, doi:10.1523/JNEUROSCI.2908-03.2004 (2004).
% matlab version R2013b

%implements the recursive filtering algorithm

%on the spike train data I



T  = size(I,2);





%sig defines size of variance in normal distn of qtk-qtk-1s

%so large sig means less constraint on continuity



qhat(1)   = qguess;    

sigsq(1)  = sige^2;   





%loop through all time

%calls newtonsolve to find solution to nonlinear

%posterior prediction estimate



%number_fail saves the time steps when Newton method fails

number_fail = [];



for t=2:T

   qhatold(t)  = nuone*qhat(t-1) ;

   sigsqold(t) = nuone*nuone*sigsq(t-1) + sige^2;



   [qhat(t), flagfail]  = newtonsolve(muone, delk, qhatold(t), sigsqold(t), I(t));

      	

   if flagfail>0

      number_fail = [number_fail t];

   end



   denom       = -1/sigsqold(t) - delk*exp(muone)*exp(qhat(t))/(1+delk*exp(muone)*exp(qhat(t)))^2 ;

   sigsq(t)    = -1/denom;

end



if isempty(number_fail)<1

   fprintf(2,'!!!!Newton convergence failed at times %d \n', number_fail)

end



p = delk*exp(muone)*exp(qhat)./(1+delk*exp(muone)*exp(qhat));





