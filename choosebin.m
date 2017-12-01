function [tmp,ind] = choosebin(logdensity,densitybins)

[tmp,ind]=min(abs(logdensity-densitybins));

end

