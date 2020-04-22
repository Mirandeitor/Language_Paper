function [timestamp,lx,ly,rx,ry] = importPopOutData(file,res,missingx,missingy,FaceChairStim)
% Imports data from Tobii TX300 as returned by Tobii SDK
% res = [xres yres]

%This function has been modified to read the excel data from the Pop-out
%task by David Lopez 28.09.2017

% FaceChairStim contains the stimuli to be analysed from the 5 faces or 5
% chairs.


[~, TXT, RAW]=xlsread(file);
%The raw data is located in column X
%The fixation data is in column Y
%Find the column with all the participant names
%First we determine the location of the columns that we are interested in
for i=1:size(RAW,2)
    switch RAW{1,i}
        case 'ParticipantName'
            participantName = i;
        case 'GazePointLeftX (ADCSpx)'
            rawFixationLeftX = i;
        case 'GazePointLeftY (ADCSpx)'
            rawFixationLeftY = i;
        case 'GazePointRightX (ADCSpx)'
            rawFixationRightX = i;
        case 'GazePointRightY (ADCSpx)'
            rawFixationRightY = i;
        case 'ValidityLeft'
            validityLeft = i;
        case 'ValidityRight'
            validityRight = i;
        case 'RecordingTimestamp'
            timestmp = i;
        case 'MediaName'
            mediaName = i;       
    end
end

stimuli.name = RAW(2,participantName);     
listMedia = unique(TXT(2:end,mediaName));      
for iMedia=1:size(listMedia,1)
    if (isempty(listMedia{iMedia}) || ~startswith(listMedia{iMedia},'Pop'))
        listMedia{iMedia} = [];
    end
end
X = RAW(2:end,mediaName);
listMedia = listMedia(~cellfun('isempty',listMedia));
for iMedia = 1:size(listMedia,1)       
	index = false(1, numel(X));
    for k = 1:numel(X)
        index(k) = strcmp(X{k,1},listMedia{iMedia,1});
    end
    stimuli.listStimuli(iMedia).name = listMedia{iMedia,1};
    diffArray = diff(index);
    stimuli.listStimuli(iMedia).endPos = find(diffArray<0);
    %It could be the case that there is not endPos since the file ends
    %in the stimuli. We have to correct for that.
    if isempty(stimuli.listStimuli(iMedia).endPos)
        stimuli.listStimuli(iMedia).endPos = size(RAW,1);
    end
    stimuli.listStimuli(iMedia).startPos = find(diffArray==1)+1;
    clear  diffArray;
end

%search for the stimuli we want and use that fixation data 
%Extract validity of right and left eye as well as fixation coordinates of
%both eyes + timestamp.
for iStimuli = 1:size(stimuli.listStimuli,2)  
    if strcmp(stimuli.listStimuli(iStimuli).name,FaceChairStim)
        valueStim = iStimuli;
        break;
    end
end

for iData = stimuli.listStimuli(valueStim).startPos:stimuli.listStimuli(valueStim).endPos
    timestamp(iData+1-stimuli.listStimuli(valueStim).startPos) = RAW{iData,timestmp};
    if isempty(RAW{iData,rawFixationLeftX})
        lx(iData+1-stimuli.listStimuli(valueStim).startPos) = missingx;
    else
        lx(iData+1-stimuli.listStimuli(valueStim).startPos) = RAW{iData,rawFixationLeftX};
    end
    if isempty(RAW{iData,rawFixationLeftY})
        ly(iData+1-stimuli.listStimuli(valueStim).startPos) = missingy;
    else
        ly(iData+1-stimuli.listStimuli(valueStim).startPos) = RAW{iData,rawFixationLeftY};
    end
    if isempty(RAW{iData,validityLeft})
        lv(iData+1-stimuli.listStimuli(valueStim).startPos) = 4;
    else
        lv(iData+1-stimuli.listStimuli(valueStim).startPos) = RAW{iData,validityLeft};
    end
    if isempty(RAW{iData,rawFixationRightX})
        rx(iData+1-stimuli.listStimuli(valueStim).startPos) = missingx;
    else
        rx(iData+1-stimuli.listStimuli(valueStim).startPos) = RAW{iData,rawFixationRightX};
    end
    if isempty(RAW{iData,rawFixationRightY})
        ry(iData+1-stimuli.listStimuli(valueStim).startPos) = missingy;
    else
        ry(iData+1-stimuli.listStimuli(valueStim).startPos) = RAW{iData,rawFixationRightY};
    end
    if isempty(RAW{iData,validityRight})
        rv(iData+1-stimuli.listStimuli(valueStim).startPos) = 4;
    else
        rv(iData+1-stimuli.listStimuli(valueStim).startPos) = RAW{iData,validityRight};
    end
end

% sometimes we have weird peaks where one sample is (very) far outside the
% monitor. Here, count as missing any data that is more than one monitor
% distance outside the monitor.
qMiss = lx<-res(1) | lx>2*res(1) | ly<-res(2) | ly>2*res(2) | lv==4;
lx(qMiss) = missingx;
ly(qMiss) = missingy;
qMiss = rx<-res(1) | rx>2*res(1) | ry<-res(2) | ry>2*res(2) | rv==4;
rx(qMiss) = missingx;
ry(qMiss) = missingy;

timestamp = timestamp';
lx = lx';
ly = ly';
rx = rx';
ry = ry';
lv = lv';
rv = rv';

return
