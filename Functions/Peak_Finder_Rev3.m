% Chooses relevant peaks out of all those returned by MATLAB. Criteria
% explained in documentation

function [Locations] = Peak_Finder_Rev3(pks,locns,p,threshold)
 
avgsnr = mean(pks);
[p,idx] = sort(p,"descend");
locns = locns(idx);
pks = pks(idx);
for i=1:length(p)
    if(p(i)<threshold)
        break
    end
end
locns = locns(1:i);
pks = pks(1:i);

[~,idx] = sort(pks,"descend");
pks=pks(idx);
locns = locns(idx);

idx = find(locns>0);
locns = locns(idx);
pks = pks(idx);

idx = find(pks>avgsnr);
locns = locns(idx);
pks = pks(idx);
%locns;
%Locations = locns(:);


% [~,idx] = sort(locns,"descend");
% if(locns(idx(1))==locns(1))
%     Locations(1) = locns(1);
% else
%     Locations(1) = locns((idx(1)));
% end
% k = 2;
% for i = 2: length(locns)
%     if(locns(i)<locns(i-1))
%         Locations(k) = locns(i);
%         k = k+1;
%     end
% end

% [~,idx] = sort(locns,"descend");
% if(locns(idx(1))==locns(1))
%     Locations(1) = locns(1);
% else

if(length(locns)>0)
    Locations(1) = locns((1));
    ref_pk = pks(1);
    ref_dist = Locations(1);
    k = 2;
    for i = 2:length(locns)
        max_pk = ref_pk - 80*log10(locns(i)/ref_dist)+10;
        min_pk = ref_pk - 80*log10(locns(i)/ref_dist)-10;
        if(pks(i)<max_pk && pks(i)>min_pk)
            Locations(k) = locns(i);
            k = k+1;
        end
    end
else
    Locations = [0];
end

% Locations(1) = locns((1));
% k = 2;
% 
% for i = 2:length(locns)
%     if(locns(i)>Locations(k-1))
%         Locations(k) = locns(i);
%         k = k+1;
%     end
% end


% if (length(locns)>numpeaks)
%     Locations = locns(1:numpeaks);
% else
%     Locations = locns(:);
% end


