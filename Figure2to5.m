% the following is an alternative way to define the 10...30 stopping time example
L=1:30;
d=[.05 .05] % change here the mortality for [females males]
competingstrains=999;
densitydependence=0;% set here wheter c should be density dependent (1) or constant (0) 
geomean=1; % set wheter you want to use geometric (1) mean or arithmetic (0) mean as optimising criterion
for rep=1:10 % How many repetition you want to do

rep


pL=ones(size(L)); pL(1:29)=0; % risky times with positive ending probability starts from cycle 10 onwards
[s1,m1,rounds,Wbest1]=genalg(L,pL,d,competingstrains,densitydependence,geomean);
prob1=cumsum(pL)/sum(cumsum(pL))

s1final(rep,:)=s1;m1final(rep,:)=m1;Wbest1final(rep,:)=Wbest1;


pL=ones(size(L)); pL(1:19)=0; % risky times with positive ending probability starts from cycle 10 onwards
[s2,m2,rounds,Wbest2]=genalg(L,pL,d,competingstrains,densitydependence,geomean);
prob2=cumsum(pL)/sum(cumsum(pL))

s2final(rep,:)=s2;m2final(rep,:)=m2;Wbest2final(rep,:)=Wbest2;

pL=ones(size(L)); pL(1:9)=0; % risky times with positive ending probability starts from cycle 10 onwards
[s3,m3,rounds,Wbest3]=genalg(L,pL,d,competingstrains,densitydependence,geomean);
prob3=cumsum(pL)/sum(cumsum(pL))

s3final(rep,:)=s3;m3final(rep,:)=m3;Wbest3final(rep,:)=Wbest3;

end

% find the best strategy and plot it
[val ind1]=max(Wbest1final)
[val ind2]=max(Wbest2final)
[val ind3]=max(Wbest3final)

s1max=s1final(ind1,:);m1max=m1final(ind1,:);
s2max=s2final(ind2,:);m2max=m2final(ind2,:);
s3max=s3final(ind3,:);m3max=m3final(ind3,:);



if densitydependence == 1

    
s=[s1final(ind1,:); s1final(ind1,:)];m=[m1final(ind1,:); m1final(ind1,:)]
[F1,M1,st1max,mt1max]=dynamics(s,m,d,[1; competingstrains],densitydependence,30);

s=[s2final(ind1,:); s2final(ind1,:)];m=[m2final(ind1,:); m2final(ind1,:)]
[F2,M2,st2max,mt2max]=dynamics(s,m,d,[1; competingstrains],densitydependence,30);

s=[s3final(ind1,:); s3final(ind1,:)];m=[m3final(ind1,:); m3final(ind1,:)]
[F3,M3,st3max,mt3max]=dynamics(s,m,d,[1; competingstrains],densitydependence,30);

figure(1); 
subplot(3,1,1); yyaxis left;bar(L,[st1max(1,:); (1-st1max(1,:)).*mt1max(1,:)]','stacked');yyaxis right; plot(1:30,sum(M1+F1));hold on; yyaxis left; plot(cumsum(prob1));hold off;
subplot(3,1,2); yyaxis left; bar(L,[st2max(1,:); (1-st2max(1,:)).*mt2max(1,:)]','stacked');yyaxis right; plot(1:30,sum(M2+F2));hold on; yyaxis left; plot(cumsum(prob2));hold off;
subplot(3,1,3); yyaxis left; bar(L,[st3max(1,:); (1-st3max(1,:)).*mt3max(1,:)]','stacked');yyaxis right; plot(1:30,sum(M3+F3));hold on; yyaxis left; plot(cumsum(prob3));hold off;


       
else

s1=[s1max;s1max];m1=[m1max;m1max];
s2=[s2max;s2max];m2=[m2max;m2max];
s3=[s3max;s3max];m3=[m3max;m3max];
[Etotal,E,F1,M1]=mutantfitness(s1,m1,d,competingstrains,densitydependence);
[Etotal,E,F2,M2]=mutantfitness(s2,m2,d,competingstrains,densitydependence);
[Etotal,E,F3,M3]=mutantfitness(s3,m3,d,competingstrains,densitydependence);

figure(1); 
subplot(3,1,1); yyaxis left;bar(L,[s1max(1,:); (1-s1max(1,:)).*m1max(1,:)]','stacked');yyaxis right; plot(1:30,sum(M1+F1));hold on; yyaxis left; plot(cumsum(prob1));hold off;
subplot(3,1,2); yyaxis left; bar(L,[s2max(1,:); (1-s2max(1,:)).*m2max(1,:)]','stacked');yyaxis right; plot(1:30,sum(M2+F2));hold on; yyaxis left; plot(cumsum(prob2));hold off;
subplot(3,1,3); yyaxis left; bar(L,[s3max(1,:); (1-s3max(1,:)).*m3max(1,:)]','stacked');yyaxis right; plot(1:30,sum(M3+F3));hold on; yyaxis left; plot(cumsum(prob3));hold off;

end








