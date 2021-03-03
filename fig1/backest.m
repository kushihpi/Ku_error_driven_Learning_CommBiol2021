function [betterq, bettersigsq, qnew, signewsq, a] = backest(q, qold, sigsq, sigsqold, nuone);
% Written by Anne C. Smith 2004
% publication: Smith, A. C. et al. Dynamic analysis of learning in behavioral experiments. The Journal of neuroscience : the official journal of the Society for Neuroscience 24, 447-461, doi:10.1523/JNEUROSCI.2908-03.2004 (2004).
% matlab version R2013b


T = size(q,2);



qnew(T)     = q(T);

signewsq(T) = sigsq(T);

for i = T-1 :-1: 1

   a(i)        = nuone*sigsq(i)/sigsqold(i+1);

   qnew(i)     = q(i) + a(i)*(qnew(i+1) - qold(i+1));

   signewsq(i) = sigsq(i) + a(i)*a(i)*(signewsq(i+1)-sigsqold(i+1));

  % signewsq(i) = sigsq(i) + a(i)*a(i)*sigsqold(i+1);

end



betterq = qnew(1);

bettersigsq = signewsq(1);

