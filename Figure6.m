
%in unpredictable environment:
L=1:30;
d=[.05 .05];
competingstrains=999;
x=[1 6 14 22 28];% Tmin
directcue=1; %choose here what cue to use; 1= direct;0=indirect

for r=1:length(x)
    r
    for rep=1:10
        [r rep]
        pL=ones(size(L)); pL(1:(30-x(r)))=0;
        
        geomean=1; % for optimisation for the geometric mean
        densitydependence=0;% change here wheter c is density depenedent(1) or not(0)
        
        if directcue == 1
            [s1,m1,rounds(r,rep,1)]=genalg(L,pL,d,competingstrains,densitydependence,geomean);
        else
            [s1,m1,rounds(r,rep,1)]=genalg(L,pL,d,competingstrains,densitydependence,geomean,30);
        end
        s=[s1; s1];m=[m1; m1];
        [Emax,allEg(rep,:)]=mutantfitness(s,m,d,competingstrains,densitydependence);
        
        
        
        [Ave(1,:,r) ind]=max(mean(allEg,2));
        Var(1,:,r)=var(allEg(ind,:));
        
        geomean=0 %Optimisation for aritmetic mean
        
        if directcue == 1
            [s1,m1,rounds(r,rep,2)]=genalg(L,pL,d,competingstrains,densitydependence,geomean);
        else
            [s1,m1,rounds(r,rep,2)]=genalg(L,pL,d,competingstrains,densitydependence,geomean,30);
        end
        s=[s1; s1];m=[m1; m1]
        
        
        [Emax,allEa(rep,:)]=mutantfitness(s,m,d,competingstrains,densitydependence);
        [Ave(2,:,r) ind]=max(mean(allEa,2));
        Var(2,:,r)=var(allEa(ind,:));
    end
   
end


%plot it

b=linspace(0,1,length(x))

for r=1:length(x)
    
    figure(2);subplot(2,2,1);plot(Ave(:,:,r),Var(:,:,r),'color',[0 b(r) 1-b(r)]); hold on
    figure(2);subplot(2,2,1);plot(Ave(1,:,r),Var(1,:,r),'ok',Ave(2,:,r),Var(2,:,r),'xk'); hold on 
  
end

