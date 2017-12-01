function  CS  = clutchsize(logdensity)% note that CS is here the log
% find clutchsize depending on current population density to implement
% density dependence in the model.

% values to define relationship from field data


CS=exp((4+(-5 -4)./(1+exp((- 0.3).*(logdensity-8)))));


%plot(logdensity,CS)


end




