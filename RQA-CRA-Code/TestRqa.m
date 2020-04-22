% This test program computes a few recurrence analyses

% Example 1:
% These are simply categories. They could, for example, be identifiers of grid elements
% on a fixation grid.

fixations=[37 28 28 28 19 18 9 1 1 9 9 1 1 9 9 17 34 42 43 44 44 45 36 28 19 11 1];

param.delay = 1;
param.embed = 1;
param.rescale = 0;
param.metric = 'euclidian';
param.adjacency=[];
param.linelength = 2;
param.radius=0.1;

result=Rqa(fixations,param);
PlotRecurrenceMatrix(result.recmat,'results/rqa1.png');

display(...
['nRec=' num2str(result.nrec) ' %Rec=' num2str(result.rec) ...
 ' Det=' num2str(result.det) ' MeanLine=' num2str(result.meanline) ...
 ' MaxLine=' num2str(result.maxline) ' ENT=' num2str(result.ent) ...
 ' RelENT=' num2str(result.relent) ' Trend=' num2str(result.trend) ...
 ' LAM=' num2str(result.lam) ' TT=' num2str(result.tt) ...
 ' corm=' num2str(result.corm)] )

pause 

% Example 2:
% The input consists of a series of fixations in (x,y) coordinates.

fixations(:,1) = [387
   190
   421
   503
   243
   185
   187
   200
   531
   496
   250
   570
   510
   234
   349
   539
   557
   302
   180
   175
   495
   507
   238
   168
   320
   488]

fixations(:,2) = [288
   348
   501
   486
   442
   312
   308
   255
   496
   476
   480
   340
   445
   285
   309
    87
   135
    77
   132
   331
   282
   493
   411
   327
   198
   158]

param.delay = 1;
param.embed = 1;
param.rescale = 0;
param.metric = 'euclidian';
param.adjacency=[];
param.linelength = 2;
param.radius=16;

result=Rqa(fixations,param);
PlotRecurrenceMatrix(result.recmat,'results/rqa2.png');

display(...
['nRec=' num2str(result.nrec) ' %Rec=' num2str(result.rec) ...
 ' Det=' num2str(result.det) ' MeanLine=' num2str(result.meanline) ...
 ' MaxLine=' num2str(result.maxline) ' ENT=' num2str(result.ent) ...
 ' RelENT=' num2str(result.relent) ' Trend=' num2str(result.trend) ...
 ' LAM=' num2str(result.lam) ' TT=' num2str(result.tt) ...
 ' corm=' num2str(result.corm)] )
