fix1=[1 2 3 4 5 6 7 8 9 10];
fix2=[4 5 6 7 8 9 10 11 12 13];

param.delay = 1;
param.embed = 1;
param.rescale = 0;
param.metric = 'euclidian';
param.adjacency=[];
param.linelength = 2;
param.radius=0.1;
param.plot=1;

result=Cra(distanceMother,distanceInfant,param);

display(...
['nRec=' num2str(result.nrec) ' %Rec=' num2str(result.rec) ...
 ' %Det=' num2str(result.det) ' MeanLine=' num2str(result.meanline) ...
 ' MaxLine=' num2str(result.maxline) ' Ent=' num2str(result.ent) ...
 ' RelEnt=' num2str(result.relent) ...
 ' Corm=' num2str(result.corm) ...
 ' Lam=' num2str(result.lam) ' TT=' num2str(result.tt)]);

PlotRecurrenceMatrix(result.recmat,'RecurrenceMatrix');
