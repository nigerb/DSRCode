function [F,M,st,mt,c,ind,binuse]=dynamics(s,m,d,F0,densitydependence,seasonlength_modelB)
% within-season dynamics (directly developing offspring only) when the season length is 
% - the length of s if seasonlength_modelB has been left undefined
% - seasonlength_modelB if this has been defined
% s = probability of going for sex
% m = probability of producing a male clutch (conditional on the cycle being asexual)
% c is clutch size (measured as those who mature)
% d(1) is female mortality, d(2) is male mortality
% F0 is the initial number of females in the population (we assume no males exist initially)

F=F0;
M=0*F;
if nargin==5 % model A has been defined
    st=s; mt=m;
    for i=1:size(s,2)
        
        if densitydependence==1
        c(i)=clutchsize(log(sum(F(:,i)+M(:,i))));
        else
            c(i)=10;
        end
        F(i+2)=F(i)*c(i)*(1-s(i))*(1-m(i));
        M(i+2)=F(i)*c(i)*(1-s(i))*m(i);
        % per step mortality is dm for males and df for females
        F(i+1)=F(i+1)+(1-d(1))*F(i);
        M(i+1)=M(i+1)+(1-d(2))*M(i);
    end
else % model B; now s, m, F0 typically have 2 rows (but also works with 1 or more than 2)
    densitybins=logspace(0,2,100);
    for i=1:seasonlength_modelB
        % first figure out current bin to use
        
        logdensity= abs(log(sum(F(:,i)+M(:,i))));
        [tmp,ind] = choosebin(logdensity,densitybins);
    
        binuse(i)=ind; 
        % then proceed as in model A
        st(:,i)=s(:,ind);
        mt(:,i)=m(:,ind);
        
        if densitydependence==1
        c(i)=clutchsize(log(sum(F(:,i)+M(:,i))));
        else
            c(i)=10;
        end
        F(:,i+2)=F(:,i).*c(i).*(1-s(:,ind)).*(1-m(:,ind));
        M(:,i+2)=F(:,i).*c(i).*(1-s(:,ind)).*m(:,ind);
        % per step mortality is dm for males and df for females
        F(:,i+1)=F(:,i+1)+(1-d(1))*F(:,i);
        M(:,i+1)=M(:,i+1)+(1-d(2))*M(:,i);
    end
end
% because the season actually ends at length(st), the beyond-that M and F
% have to be removed (they never had the time to mature in reality)
F=F(:,1:length(st));
M=M(:,1:length(st));

