function [RRNormalised, fixationRecMatrix ] = obtainNormalisedRecValues(fixationX, fixationY, duration, recMatrix,params)

% This function is going to obtain the normalised Recurrence values from
% the recurrence matrix using the fixation in X and y coordinates, in
% combination with the fixation durations.

% Input parameters: 
    % - fixationX and fixationY are the x and y coordinates of the fixations.
    % - duration are the fixation duration for each XY combination
    % - recMatrix contains the recurrences of the spatial fixations
    
% Output parameters:    
    % - RRNormalised -> recurrence values spatially normalised taking into account fixations
    % durations.
    % - fixationRecMatrix -> recurrence matrix containing the fixation
    % durations.

%V1.0 Creation of the document by David Lopez 20.10.2016
%V1.1 Bug Fix when no determinism or laminarity was foung the value should
%be 0 not NaN by David Lopez 28.10.2016
%V1.2 Bug Fix the Rec now returns 0 instead of NaN when no recurrences are
%found 28.10.2016
%V1.3 Bug Fix The CORM values was being wrong. Now the values have been
%corrected by David Lopez 28.10.2016

if nargin<4
    error('The recurrence matrix has not been provided');
end

if nargin<3
    error('The duration array has not been provided');
end

if nargin<2
    error('The fixationY array has not been provided');
end

if nargin<1
    error('The fixationX array has not been provided');
end

%default params
if nargin<5
    minLineLength
end

%Initialise variables
fixationRecMatrix = double(zeros(size(recMatrix)));

for iFixX = 1:size(fixationX,1)
    for iFixY = 1:size(fixationY,1)
        % We compare each x position with all the y in the plot. If there
        % is a recurrence we transform that position into the fixation i +
        % fixation j value.
        if recMatrix(iFixX,iFixY) == 1
            fixationRecMatrix(iFixX,iFixY) = duration(1,iFixX) + duration(1,iFixY);
        else
            continue;
        end
    end
end
%Calculate the recurrence parameters from fixationRecMatrix
overallDuration = sum(duration);
numFixations = size(fixationX,1);
%Data from the main diagonal is not needed
partialRecurrenceMatrix=triu(fixationRecMatrix,1);

%Normalised Recurrence
if sum(sum(partialRecurrenceMatrix)) == 0
    RRNormalised.RECn = 0;
else    
    RRNormalised.RECn = 100*sum(sum(partialRecurrenceMatrix))/((numFixations-1)*overallDuration);%Otherwise that can give NaN
end
% Diagonals of recurrence matrix, excluding main diagonal. 
[recurrenceDiagonals,diagonals]=spdiags(partialRecurrenceMatrix);
recurrenceDiagonalsAux = recurrenceDiagonals;
recurrenceDiagonalsAux(recurrenceDiagonalsAux>0)=1;
nDiagonals=length(diagonals);
thresholdedDiagonals=[];
for i=1:nDiagonals
    d=[0; recurrenceDiagonalsAux(:,i); 0];
    dchange=diff(d);
    dlengthNeg = find(dchange==-1);
    dlengthPos = find(dchange==1);
    for iLength = 1:(size(dlengthNeg,1))
        if (dlengthNeg(iLength,1) - dlengthPos(iLength,1)) >= params.linelength
            thresholdedDiagonals=[thresholdedDiagonals; sum(recurrenceDiagonals(dlengthPos(iLength):(dlengthNeg(iLength)-1),i))];
        else
            continue;
        end
    end
end
%Normalised Determinism
if ~isempty(thresholdedDiagonals)
   RRNormalised.DETn =100*sum(thresholdedDiagonals)/sum(sum(partialRecurrenceMatrix));  
else
   RRNormalised.DETn = 0;
end
CORMvalues = [];
for i = 1:(numFixations-1)
    for j = (i+1):numFixations
        CORMvalues = [CORMvalues;(j-i)*partialRecurrenceMatrix(i,j)];
    end
end
RRNormalised.CORMn=100*sum(CORMvalues)/(((numFixations-1)^2)*overallDuration);

%Vertical and horizontal lines
thresholdedVerticals=[];
paddedRM=[zeros(numFixations,1) recMatrix-eye(numFixations) zeros(numFixations,1)];
% In CRA, diagonal elements can be zero 
paddedRM=max(paddedRM,0);
for i=1:numFixations
    v=paddedRM(i,:);
    vchange=diff(v);
    vlengthNeg=find(vchange<0);
    vlengthPos=find(vchange>0);
    if size(vlengthPos,2)>0
        for iLength = 1:(size(vlengthPos,2))
            if (vlengthNeg(1,iLength) - vlengthPos(1,iLength)) >= params.linelength
                thresholdedVerticals=[thresholdedVerticals; sum(fixationRecMatrix(i,vlengthPos(iLength):(vlengthNeg(iLength)-1)))];
            else
                continue;
            end
        end
    end
end 
%Normalised laminarity
if ~isempty(thresholdedVerticals)
    RRNormalised.LAMn=100*sum(thresholdedVerticals)/(2*sum(sum(partialRecurrenceMatrix)));
else
    RRNormalised.LAMn = 0;
end        
