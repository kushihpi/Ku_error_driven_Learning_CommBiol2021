function [dense, inst] = f_my_smoothfr (spmat, kwidth, ktype)
%SMOOTHFR  smoothed firing rates
% ex: [dense, inst] = f_my_smoothfr (list_firing_rate, 5, 'gauss');
% make sure the kwidth is of even size
% spmat: nx1 array
kwidth = kwidth - rem (kwidth, 2);
halfwin = kwidth/2;
[inst, spl] = size (spmat);

if  inst == 0, dense=[]; return, end;
%if inst ==1 then average firing rates assumed
if  inst >  1, spmat = sum (spmat, 1)/inst*1000; end;

switch lower(ktype)
case 'gauss'
   x = -halfwin:halfwin;
   kernel = exp (-x.^2/halfwin^2);
case 'boxcar'
   kernel = ones (kwidth + 1, 1);
case 'epsp'
   t=0:kwidth;
   tg=1;
   td=20;
   r = (1-exp (-t/tg)).*(exp(-t/td));
   kernel = r /sum(r);
otherwise
   error ('type should be one of ''gauss'', ''boxcar'', ''epsp''');
end
kernel = kernel / sum (kernel);
padded = zeros (1, spl + kwidth * 2);
padded (halfwin + 1 : spl + halfwin) = spmat;
padded (1:halfwin) = repmat (mean (spmat (1:halfwin)), 1, halfwin);
padded (halfwin + spl + 1 : spl + kwidth) = repmat (mean (spmat(spl-halfwin:spl)), 1, halfwin);
dense = filter (kernel, 1, padded);
dense = dense (kwidth+1 : end-kwidth);
% figure
% plot(dense)
