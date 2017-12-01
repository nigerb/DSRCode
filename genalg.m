function [s,m,rounds,Wbest,st,mt,binuse]=genalg(L,pL,d,competingstrains,densitydependence,geomean,modelB)

% L contains the season lengths that are possible
% pL indicates the corresponding probability that the season ends at L


p=0.01; % this is the probability of copying a neighbouring strategy

% initial strategies are uniformly distributed random numbers between 0 and 1
if nargin==6 % model A
    s_resident=rand([1 max(L)]);
    m_resident=rand([1 max(L)]);
else % model B is always run with 30 density bins
    s_resident=rand([1 100]);
    m_resident=rand([1 100]);
  
end

% let the genetic algorithm run for .... generations
for t=1:80000
     s=[s_resident; s_resident]; m=[m_resident; m_resident];
     if nargin== 6
       [Etotal,E,F,M,c]=mutantfitness(s,m,d,competingstrains,densitydependence);
       
        

     else
     [F,M,st,mt,c,ind,binuse]=dynamics(s,m,d,[1; competingstrains],densitydependence,modelB);
     %%%

     end
     
  
     
    % we test the 99 novel strategies against the resident and store their fitness in W(1), W(2), ..., W(9)
    Smut=ones([99 1])*s_resident+0.001*randn([99 length(s_resident)]);
    Mmut=ones([99 1])*m_resident+0.001*randn([99 length(m_resident)]);
    Smut(Smut<0)=0; Smut(Smut>1)=1;
    Mmut(Mmut<0)=0; Mmut(Mmut>1)=1;

    if nargin==6
        for i=1:99
            % each variant is tested against the resident
            s=[Smut(i,:); s_resident];
            m=[Mmut(i,:); m_resident];
            
            if geomean==1
            W(i)=loggeommean(L,pL,s,m,d,competingstrains,densitydependence);
            else
                W(i)=aritmean(L,pL,s,m,d,competingstrains,densitydependence);
            end
        end
    else % model B
        for i=1:99
            % each variant is tested against the resident
            s=[Smut(i,:); s_resident];
            m=[Mmut(i,:); m_resident];
            
            if geomean==1
            W(i)=loggeommean(L,pL,s,m,d,competingstrains,densitydependence,1);
            else
                W(i)=aritmean(L,pL,s,m,d,competingstrains,densitydependence,1);
            end
            
        end

    end
    
    % we also need to compute the fitness of the resident against itself
    s=[s_resident; s_resident]; m=[m_resident; m_resident];
    if nargin==6
        if geomean==1
        W(100)=loggeommean(L,pL,s,m,d,competingstrains,densitydependence);
        else
                W(100)=aritmean(L,pL,s,m,d,competingstrains,densitydependence);
            end
    else
        if geomean==1
        W(100)=loggeommean(L,pL,s,m,d,competingstrains,densitydependence,1);
        else
                W(100)=aritmean(L,pL,s,m,d,competingstrains,densitydependence,1);
            end
    end
    
    % Check which variants won
    [tmp,order]=sort(-W); 
    intonextgen=[1 order(1:2)]'; % the resident and the 9 best perfoming ones form recombinants
    Schoices=[Smut; s_resident];
    Mchoices=[Mmut; m_resident];
    Srecombinant=[]; Mrecombinant=[];
    for i=1:length(s_resident)
        Srecombinant(:,i)=Schoices(intonextgen(unidrnd(3,9,1)),i);
        Mrecombinant(:,i)=Mchoices(intonextgen(unidrnd(3,9,1)),i);
    end
    
    %now occasionally copy the neighbouring strategy (with probability p)
    %...but the first generation can only copy the future strategy (probability p/2)
    whatwehavenowS=Srecombinant; % a temporary storage for strategies right now
    whatwehavenowM=Mrecombinant;
    for j=1:size(Srecombinant,1)
       if rand<p/2 Srecombinant(j,1)=mean(whatwehavenowS(j,1:2)); end; % the 1st one can only copy from the future
        for i=2:max(L)-1
            if rand<p % generations 2 to max(L) can copy from the past or the future
                if rand<0.5 Srecombinant(j,i)=mean(whatwehavenowS(j,i-1:i)); else Srecombinant(j,i)=mean(whatwehavenowS(j,i:i+1)); end
            end
        end
        if rand<p/2 Srecombinant(j,max(L))=mean(whatwehavenowS(j,max(L)-1:max(L))); end; % last generation can only copy from the past
    end
    
    for j=1:size(Mrecombinant,1)
        if rand<p/2 Mrecombinant(j,1)=mean(whatwehavenowM(j,1:2)); end; % the 1st one can only copy from the future
        for i=2:max(L)-1
            if rand<p % generations 2 to max(L) can copy from the past or the future
                if rand<0.5 Mrecombinant(j,i)=mean(whatwehavenowM(j,i-1:i)); else Mrecombinant(j,i)=mean(whatwehavenowM(j,i:i+1)); end
            end
        end
        if rand<p/2 Mrecombinant(j,max(L))=mean(whatwehavenowM(j,max(L)-1:max(L))); end; % last generation can only copy from the past
    end
    
    % Test the 9 new recombinants against the resident
    
    W(10)=W(100);
    if nargin==6
        for i=1:9
            s=[Srecombinant(i,:); s_resident];
            m=[Mrecombinant(i,:); m_resident];
            
            if geomean==1
            W(i)=loggeommean(L,pL,s,m,d,competingstrains,densitydependence);
            else
                W(i)=aritmean(L,pL,s,m,d,competingstrains,densitydependence);
            end
        end
    else
        for i=1:9
            s=[Srecombinant(i,:); s_resident];
            m=[Mrecombinant(i,:); m_resident];
            if geomean==1
            W(i)=loggeommean(L,pL,s,m,d,competingstrains,densitydependence,1);
            else
                W(i)=aritmean(L,pL,s,m,d,competingstrains,densitydependence,1);
            end
        end
    end
    
    

    % now choose the best one to be the new resident
    [Wbest,ind]=max(W);
    fit(t)=Wbest;
    
    if t>200
    if range(fit(t-200:t))==0.0000001
        rounds=t
        disp('Yay')
        break
    end
    end
    
   
    if ind<10 s_resident=Srecombinant(ind,:); m_resident=Mrecombinant(ind,:); end % if not, then we keep the resident
    
    
    rounds=t;
end
s=s_resident; m=m_resident;
 