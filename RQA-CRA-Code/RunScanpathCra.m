% Compute scanpath cra

% General parameters

inputFilename='data/cra-fixations.txt';

% Output parameters
scanpathParam.plotFixations=1;
scanpathParam.fixationRadius=64;
scanpathParam.plotRecurrenceMatrix=1;
scanpathParam.screenx=1024;
scanpathParam.screeny=768;
scanpathParam.outlierTreatment='eliminate'; % none,eliminate,include

% Cra parameters
param.delay=1;
param.embed=1;
param.rescale=0;
param.metric='euclidian';
param.linelength=2;
param.radius=64;
param.adjacency=0;

% ---------
% Read data
% ---------

% Read file and extract relevant columns from spreadsheet
inputFile=fopen(inputFilename,'r');
% Skip header
format='%s %s %s %s %s %s %s %s %s';
data=textscan(inputFile,format,1,'delimiter','\t');
% Read data
format='%d %d %d %f %f %d %d %d %s';
data=textscan(inputFile,format,'delimiter','\t','treatAsEmpty','.');
fclose(inputFile);

dataTrialNumber=data{1};
dataFirstSecond=data{2};
dataFixationX=data{4};
dataFixationY=data{5};
dataDuration=data{6};
dataAOI=data{7};
dataImageNumber=data{8};
dataImageName=data{9};

% Read data
nImages=length(unique(dataImageNumber));
nLines=size(dataTrialNumber,1);
data=cell(nImages,2);

startTrial=dataTrialNumber(1);
for iLine=2:nLines+1

    if iLine<=nLines && dataTrialNumber(iLine)==dataTrialNumber(iLine-1)
        continue
    end

    endTrial=iLine-1;

    firstSecond=dataFirstSecond(startTrial);
    nFixations=endTrial-startTrial+1;
    x=dataFixationX(startTrial:startTrial+nFixations-1);
    y=dataFixationY(startTrial:startTrial+nFixations-1);
    fixations=[x y];
    imageNumber=dataImageNumber(startTrial);

    data{imageNumber,firstSecond}.fixations=fixations;
    data{imageNumber,firstSecond}.imageName=dataImageName{startTrial};

    startTrial=iLine;
end 

% Analyze data

for imageNumber=1:nImages

    fixation1=data{imageNumber,1}.fixations;
    fixation2=data{imageNumber,2}.fixations;
    imageName=data{imageNumber,1}.imageName;

    % Trim data to equal length
    n1=size(fixation1,1);
    n2=size(fixation2,1);
    if n1 ~= n2
        n=min(n1,n2);
        fixation1=fixation1(1:n,:);
        fixation2=fixation2(1:n,:);
    end

    % Plot fixations
    if scanpathParam.plotFixations
        imageFilename=['images/' imageName];
        PlotFixationsCraDistance(fixation1,fixation2,64,imageFilename);
    end

    % Compute Cra
    result = Cra(fixation1,fixation2,param);

    display(...
    ['nRec=' num2str(result.nrec) ' %Rec=' num2str(result.rec) ...
     ' %Det=' num2str(result.det) ' MeanLine=' num2str(result.meanline) ...
     ' MaxLine=' num2str(result.maxline) ' Ent=' num2str(result.ent) ...
     ' RelEnt=' num2str(result.relent) ...
     ' Corm=' num2str(result.corm) ...
     ' Lam=' num2str(result.lam) ' TT=' num2str(result.tt)]);

    if scanpathParam.plotRecurrenceMatrix
        plotFilename=['results/cra-' imageName];
        PlotRecurrenceMatrix(result.recmat,plotFilename);
    end

    pause
end
