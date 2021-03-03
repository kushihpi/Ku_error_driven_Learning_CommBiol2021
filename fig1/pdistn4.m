function [p05, p95, pmid, pmode] = pdistn4(q, s, muone, delk);
% Written by Anne C. Smith 2004
% publication: Smith, A. C. et al. Dynamic analysis of learning in behavioral experiments. The Journal of neuroscience : the official journal of the Society for Neuroscience 24, 447-461, doi:10.1523/JNEUROSCI.2908-03.2004 (2004).
% matlab version R2013b

%loop through each trial
for ov = 1:size(q,2)

 qq = q(ov);
 ss = s(ov);

 dels=1e-4;% can be changed to 5 if there is too many 1's in a row. but it does make the program run much slower. 
 

 %computes pdf at points between 0 and 1
 pr  = dels:dels:1-dels;

 %uses change of variable formula to get distn of probability

 fac   = log( pr./(1-pr)/exp(muone)/delk ) - qq;
 fac   = exp(-fac.^2/2/ss);
 pd    = dels*(sqrt(1/2/pi/ss) * 1./(pr.*(1-pr)).* fac);

 sumpd = cumsum(pd);

 %find 95 percent CLs
% lowlimit  = find(sumpd>0.025);
 lowlimit  = find(sumpd>0.05);

 lowlimit  = lowlimit(1);
% highlimit = find(sumpd>0.975);
 highlimit = find(sumpd>0.95);

 highlimit = highlimit(1)-1;
 middlimit = find(sumpd>0.5);
 middlimit = middlimit(1);

 p05(ov)   = pr(lowlimit(1));
 p95(ov)   = pr(highlimit(1)-1);
 pmid(ov)  = pr(middlimit(1));
 [y,i]     = max(pd);
 pmode(ov) = pr(i);


end


