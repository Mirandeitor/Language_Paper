function [data] = extractDataforI2MC(data,processingData,res,missingx,missingy,freq,codes)
%This function is going to return the neccesary data
%to extract fixations using the I2MC algorithm


%Validation of the input parameters
if  ((nargin<7) || (isempty(codes)))
    warning('The column names are needed');
end


if ((nargin<6) || (isempty(freq)))
    ETfreq = 1/120;%T120 assumed
else
    ETfreq = 1/freq;
end

if ((nargin<5) || (isempty(missingy)))
    missingy = -1;
end

if ((nargin<4) || (isempty(missingx)))
    missingx = -1;
end

if ((nargin<3) || (isempty(freq)))
    res = [1280 1024];%T120 assumed
end

if ((nargin<2) || (isempty(processingData)))
    error('No processed data has been included');
end

if ((nargin<1) || (isempty(data)))
    error('No  data has been included');
end

%search for the stimuli we want and use that fixation data 
%Extract validity of right and left eye as well as fixation coordinates of
%both eyes + timestamp

for iStimuli= 1:size(data.stimuli,2)
    for iTrial = 1:size(data.stimuli(iStimuli).reactionTimes,2)
        for iData = data.stimuli(iStimuli).startPos(1,iTrial):data.stimuli(iStimuli).endPos(1,iTrial)
            %Timestampt Info
            try
                timestamp(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = str2num(processingData.textdata{iData,codes.timestamp});
            catch
                timestamp(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = -1;
            end
            %Raw x left eye
            try
                if isempty(processingData.textdata{iData,codes.GazePointLeftX})
                    lx(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = missingx;
                else                
                    lx(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = str2num(processingData.textdata{iData,codes.GazePointLeftX});
                end
            catch
                lx(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = missingx;
            end
       
            %Raw y left eye
            try
                if isempty(processingData.textdata{iData,codes.GazePointLeftY})
                    ly(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = missingy;
                else
                    ly(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = str2num(processingData.textdata{iData,codes.GazePointLeftY});
                end
            catch
                ly(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = missingy;
            end
            
            %Validation left eye
            try
                if (isempty(processingData.data(iData,codes.ValidityLeft)) || isnan(processingData.data(iData,codes.ValidityLeft)))
                    lv(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = 4;
                else
                    lv(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = processingData.data(iData-1,codes.ValidityLeft);
                end
            catch
                lv(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = 4;
            end
            %Raw x right eye
            try
                if isempty(processingData.textdata{iData,codes.GazePointRightX})
                    rx(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = missingx;
                else
                    rx(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = str2num(processingData.textdata{iData,codes.GazePointRightX});
                end
            catch
                rx(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = missingx;
            end
            %Raw y right eye
            try
                if isempty(processingData.textdata{iData,codes.GazePointRightY})
                    ry(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = missingy;
                else                
                    ry(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = str2num(processingData.textdata{iData,codes.GazePointRightY});
                end
            catch
            	ry(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = missingy;
            end
            %Validation right eye
            try
                if (isempty(processingData.data(iData,codes.ValidityRight)) || isnan(processingData.data(iData,codes.ValidityRight)))
                    rv(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = 4;
                else            
                    rv(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = processingData.data(iData-1,codes.ValidityRight);
                end
            catch
                rv(iData+1-data.stimuli(iStimuli).startPos(1,iTrial)) = 4;
            end                        
        end
        % Sometimes we have weird peaks where one sample is (very) far outside the
        % monitor. Here, count as missing any data that is more than one monitor
        % distance outside the monitor.
        qMiss = lx<-res(1) | lx>2*res(1) | ly<-res(2) | ly>2*res(2) | lv==4;
        lx(qMiss) = missingx;
        ly(qMiss) = missingy;
        qMiss = rx<-res(1) | rx>2*res(1) | ry<-res(2) | ry>2*res(2) | rv==4;
        rx(qMiss) = missingx;
        ry(qMiss) = missingy;
        data.stimuli(iStimuli).lx{iTrial} = lx';
        data.stimuli(iStimuli).ly{iTrial} = ly';
        data.stimuli(iStimuli).rx{iTrial} = rx';
        data.stimuli(iStimuli).ry{iTrial} = ry';
        data.stimuli(iStimuli).lv{iTrial} = lv';
        data.stimuli(iStimuli).rv{iTrial} = lv';
        data.stimuli(iStimuli).timestamp{iTrial} = timestamp';
        %Clear variables
        clear timestamp lx ly lv rx ry rv
    end
    
end
