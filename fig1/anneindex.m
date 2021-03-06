%function [learning_trial,middle,sm05, sm95]=anneindex(list)  
% Written by Anne C. Smith 2004
% publication: Smith, A. C. et al. Dynamic analysis of learning in behavioral experiments. The Journal of neuroscience : the official journal of the Society for Neuroscience 24, 447-461, doi:10.1523/JNEUROSCI.2908-03.2004 (2004).
% matlab version R2013b
list=[0 1 0 1 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1]

 I=list ; %important list HAS to be horizontal
 replace=I-3;
 replace=replace';

check_sum=sum(I);
        total=length(I);
        
        delk        = 1;
        
        sige        = 0.6;       %this value may be varied
        
        rhoone      = 0;
        
        qguess      = 0;
        
        nuone       = 1.00;   
        
        %muone       = -log(3);   %sets background probability to 0.25 
        %muone=-log(2); %set proba to 0.33
         muone = 0.5; %sets background probability to 0.5 
        
        %call recfilter to filter data and return estimates for q and variance
                
        [p, q, s, qold, sold]=recfilter(I, sige, qguess,rhoone, delk, nuone, muone);
        
        [betterq, bettersigsq, qnew, signewsq, a]=backest(q, qold, s, sold, nuone);
              
        %call pdistn (uses change of variable formula) to estimate conf limits of probability
        
        try
            %  [p05, p95, pmid, pmode] = pdistn4(q, s, muone, delk);
            [sm05, sm95, smmid, smmode] = pdistn4(qnew, signewsq, muone, delk);
            
        catch
            % [p05, p95, pmid, pmode] = pdistn5(q, s, muone, delk);
           try
               [sm05, sm95, smmid, smmode] = pdistn5(qnew, signewsq, muone, delk);
           catch 
               [smmid]=[replace]; [sm05]=replace;
           end
        end
        % [sm05, sm95, smmid, smmode] = pdistn(qnew, signewsq, muone, delk);
        
         middle=[smmid];       
        lc = [sm05];
        
        
        h=lc<0.33;
        if (sum(h)==length(lc))
            index=-1;
        else
            index=sum(h)+1;
        end
        learning_trial=index
  
learning_curves=[sm05; sm95; smmid]'
plot(learning_curves)