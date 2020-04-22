% Compute cra of fixation data

% General parameters
inputFilename='data/cra-fixations.txt';
outputFilename='results/cra-direct-results.txt';
outlierTreatment='none'; % none,eliminate,include

% Fixation plots
plotFixations=0;
fixationRadius=64;

% Cra parameters
plotRqaMatrix=0;
param.delay=1;
param.embed=1;
param.rescale=0;
param.metric='euclidian';
param.linelength=2;
param.radius=64;
param.adjacency=0;

% screen parameters
screenx=1024;
screeny=768;

% Open output file and write header
output=fopen(outputFilename,'w');
fprintf(output,'trial\timage\trec\tdet\tlam\ttt\tcorm\tent\trelent\n');
fprintf(1,'trial\timage\trec\tdet\tlam\ttt\tcorm\tent\trelent\n');

% Read file and extract columns
inputFile=fopen(inputFilename,'r');
% Skip header
textscan(inputFile,'%s %s %s %s %s',1);
% Read rest of data file
data=textscan(inputFile,'%d %s %f %f %d');
fclose(inputFile);

% Extract data columns
trialNumber=data{1};
imageName=data{2};
fixationX=data{3};
fixationY=data{4};
fixationDur=data{5};
nLines=size(trialNumber,1);   

startTrial=1;
for iLine=2:nLines+1

    if iLine<=nLines && trialNumber(iLine)==trialNumber(iLine-1)
        continue
    end
    
    endTrial=iLine-1;
    trial=trialNumber(startTrial);
    x=fixationX(startTrial:endTrial);
    y=fixationY(startTrial:endTrial);
    fixation=[x y];
    name=imageName{startTrial};

    % Plot fixations
    if plotFixations==1
        fixationPlotFilename=['results/rqa-direct-fixations-' num2str(trial)];
        imageFilename=['images/' name];
        PlotFixationsDistance(x,y,fixationRadius,imageFilename,fixationPlotFilename);
    end
    
    if ~strcmp(outlierTreatment,'none')
        outliers = find(x<0 | x>=screenx | y<0 | y>=screeny);
        if strcmp(outlierTreatment,'eliminate')
            fixation(outliers,:)=[];
        elseif strcmp(outlierTreatment,'include')  
            % Outliers are set to negative values with differences
            % larger than reasonable radius values to ensure that
            % there are no outlier recurrences
            nOutliers = length(outliers);
            minOutlierDistance=max(1000,10*param.radius);
            outlierValues=(-1:-1:-nOutliers)'*minOutlierDistance;
            fixation(outliers,:)=[outlierValues outlierValues];
        end
    end
    
    % Rqa analysis
    result=Rqa(fixation,param);
    if ~isstruct(result) || result.nrec==0
         continue
    end

    % Plot recurrence matrix
    if plotRqaMatrix==1
        recurrencePlotFilename=['results/rqa-direct-recurrence-' num2str(trial)];
        PlotRecurrenceMatrix(result.recmat,recurrencePlotFilename);
    end
    
    % Print results.
    [outputLinePart1,~]=sprintf('%d\t%s',...
        trial,name);
    [outputLinePart2,~]=sprintf('\t%.2f',...
        result.rec,result.det,result.lam,result.tt,result.corm,...
        result.ent,result.relent);
    outputLine=[outputLinePart1 outputLinePart2];
    % Replace NaN's by period
    outputLine=strrep(outputLine,'NaN','.');
    fprintf(output,'%s\n',outputLine);
    fprintf(1,'%s\n',outputLine);
    
    pause
    
    startTrial=iLine;
end
fclose(output);
