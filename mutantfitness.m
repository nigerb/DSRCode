function [Etotal,E,F,M,c]=mutantfitness(s,m,d,competingstrains,densitydependence,seasonlength_modelB);

% notation: 1st row of s and m = what the mutant will do; 
% 2nd row = what the resident population will do
% c is clutch size (measured as those who mature)
% d(1) is female mortality, d(2) is male mortality
% competingstrains gives the relative abundance of the resident population,
% e.g. if this is 999 then the total pop starts at 1000 of which one is the mutant herself.

if nargin==5
% model A makes the assumption that the horizontal length of s is the
% season length. 
    
    % mutant population dynamics, then population dynamics
    
    [Fmut,Mmut,st,mt,c]=dynamics(s(1,:),m(1,:),d,1,densitydependence);
    [Fpop,Mpop,st,mt,c]=dynamics(s(2,:),m(2,:),d,competingstrains,densitydependence);
    F(1,:)=Fmut;M(1,:)=Mmut;
    F(2,:)=Fpop;M(2,:)=Mpop;
    

else % model B
    % unlike in A, both dynamics are now computed at the same time
    % hence also add the mutant explicitly to competingstrains
    [F,M,st,mt]=dynamics(s,m,d,[1; competingstrains],densitydependence,seasonlength_modelB);
    Fmut=F(1,:); Mmut=M(1,:);
    Fpop=F(2,:); Mpop=M(2,:);
end

% now the production of eggs that go dormant (inside ephippia). 

% below, 'own' refers to mutant, 'other' to pop (because this file is meant
% to only compute fitness for mutants)
sexpossible=(Mmut+Mpop)>0; % if sexpossible(i) is 1, then eggs can be fertilized, otherwise not
if nargin==5
    Ephippia_own=Fmut.*s(1,:).*sexpossible;
    Ephippia_others=Fpop.*s(2,:).*sexpossible;
else
	Ephippia_own=Fmut.*st(1,:).*sexpossible;
    Ephippia_others=Fpop.*st(2,:).*sexpossible;
end

ownmales=Mmut./Mpop; ownmales(isnan(ownmales))=0; 
% now the gene copies present in own and others' ephippia (others can contain own genes because own males can mate with them)
E=Ephippia_own.*(1+ownmales)+Ephippia_others.*ownmales;

Etotal=sum(E);