function w=aritmean(L,pL,s,m,d,competingstrains,densitydependence,modelB)
% here L is all possible season lengths (e.g. 10:30, or 1:30, whatever you like but it should be integers)
% pL, which must be the same length of L, should be constructed like this: pL(i) denotes the probability that the season
% ends at each of the possible times L(i). Note that the function below
% makes sure that the values in pL sum up to 1.
% s, m, c and d should be constructed like in myownfitness.m; 
% note that the horizontal length of s must equal the largest (final) number of L, so
% that it is defined what individuals do at all the possible seasons.

% normalize pL
pL=pL./sum(pL);

if nargin==7
    [Emax,allE]=mutantfitness(s,m,d,competingstrains,densitydependence); % we make use of the fact that while Emax would happen if the season was as long as possible, allE gives all the values up to that point
else % model B
    [Emax,allE]=mutantfitness(s,m,d,competingstrains,densitydependence,L(end)); 
end
    
% the next bit has to be done carefully: we need to ignore those things
% that never happen (pL=0) but would lead to zero ephippia, which is done by
% realizing that those cases lead to fitness NaN, so taking the nansum
% will help
E_so_far=cumsum(allE(L));
fitnesscontributions=pL.*E_so_far;
w=nanmean(fitnesscontributions);




