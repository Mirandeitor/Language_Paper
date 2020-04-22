%this script is going to analyse the data extracted from I2MC algorithm
clear all;
%Load the data from fixations
file = '/Users/lopezd/Documents/Babylab/Files/London Eye Tracking Data/ET Raw Data Talby/Output/Chair 6/allfixations.txt';
dat = importdata(file);
%Convert the data to the format i had before
for i = 2:size(dat,1)
    cellAux = strsplit(dat{i}); 
    participantArray{i-1} = cellAux{15};
end
participantNames = unique(participantArray);
for iPart = 1:size(participantNames,2)
    listParticipants(iPart).name = participantNames{1,iPart};
    counter = 1;
    for i = 2:size(dat,1)
        cellAux = strsplit(dat{i});
        if strcmp(cellAux{15},participantNames{1,iPart})
            listParticipants(iPart).fixationX(counter) = str2num(cellAux{4});
            listParticipants(iPart).fixationY(counter) = str2num(cellAux{5});
            listParticipants(iPart).durations(counter) = str2num(cellAux{3});
            counter = counter + 1;
        end
    end
end
%Clean negative fixations and removed then from the analysis
for iPart = 1:size(listParticipants,2)
    durations = listParticipants(iPart).durations;
    fixationX = listParticipants(iPart).fixationX;
    fixationY = listParticipants(iPart).fixationY;
    arrayToDelete = (fixationX < 0) | (fixationY < 0 | fixationX > 800 | fixationY > 600);
    durations(arrayToDelete) = [];
    fixationX(arrayToDelete) = [];
    fixationY(arrayToDelete) = [];
    listParticipants(iPart).durations = durations;
    listParticipants(iPart).fixationX = fixationX;
    listParticipants(iPart).fixationY = fixationY;
    clear arrayToDelete
end

%First the names of the columns
cellToWrite{1,1,1} = 'ParticipantName';
cellToWrite{1,1,2} = 'NumberofRecurrences';
cellToWrite{1,1,3} = '%ofRecurrence';
cellToWrite{1,1,4} = 'Determinism';
cellToWrite{1,1,5} = 'MeanLine';
cellToWrite{1,1,6} = 'MaxLine';
% cellToWrite{1,1,7} = 'Entropy';
cellToWrite{1,1,8} = 'RelEnt.';
cellToWrite{1,1,9} = 'Trend';
cellToWrite{1,1,10} = 'Laminarity';
cellToWrite{1,1,11} = 'TrappingTime';
cellToWrite{1,1,12} = 'Corm';
cellToWrite{1,1,13} = 'RRNormalised';
cellToWrite{1,1,14} = 'DETNormalised';
cellToWrite{1,1,15} = 'CORMNormalised';
cellToWrite{1,1,16} = 'LAMNormalised';
cellToWrite{1,1,17} = 'NumberOfFixations';
%Start the RQA process
param.delay = 1;
param.embed = 1;
param.rescale = 0;
param.metric = 'euclidian';
param.adjacency=[];
param.linelength = 4;
param.radius = 64;
for participantNumber = 1:size(listParticipants,2)
	fixations(:,1) =  listParticipants(participantNumber).fixationX';
	fixations(:,2) =  listParticipants(participantNumber).fixationY';
    %Check that there are at least 5 fixations
    if size(fixations,1)>4 
    	result=Rqa(fixations,param);
        resultFace6{participantNumber}=RqaForHeatMap(fixations,param);
        %fileSave = strcat('/Users/lopezd/Documents/Babylab/Files/London Eye Tracking Data/Tobi studio data/TALBY1-Popout_Pop out_Face 6/',participantList(participantNumber).Name);
        %PlotRecurrenceMatrix(result{participantNumber}.recmat,fileSave);
        cellToWrite{participantNumber+1,1} = listParticipants(participantNumber).name;
        cellToWrite{participantNumber+1,2} = num2str(result.nrec);
        cellToWrite{participantNumber+1,3} = num2str(result.rec);
        cellToWrite{participantNumber+1,4} = num2str(result.det);
        cellToWrite{participantNumber+1,5} = num2str(result.meanline);
        cellToWrite{participantNumber+1,6} = num2str(result.maxline);
        cellToWrite{participantNumber+1,7} = num2str(result.ent);
        cellToWrite{participantNumber+1,8} = num2str(result.relent);
        cellToWrite{participantNumber+1,9} = num2str(result.trend);
        cellToWrite{participantNumber+1,10} = num2str(result.lam);
        cellToWrite{participantNumber+1,11} = num2str(result.tt);
        cellToWrite{participantNumber+1,12} = num2str(result.corm);
            
        [RRNormalised, fixationRecMatrix ] = obtainNormalisedRecValues(listParticipants(participantNumber).fixationX',listParticipants(participantNumber).fixationY', ...
                listParticipants(participantNumber).durations, result.recmat,param);
        cellToWrite{participantNumber+1,13} = num2str(RRNormalised.RECn);
        cellToWrite{participantNumber+1,14} = num2str(RRNormalised.DETn);
        cellToWrite{participantNumber+1,15} = num2str(RRNormalised.CORMn);
        cellToWrite{participantNumber+1,16} = num2str(RRNormalised.LAMn);
        cellToWrite{participantNumber+1,17} = num2str(size(fixations,1));        
            %{
            plotFilename = strcat(PATHSTR,'/',participantList(participantNumber).Name,'PLOTDISTANCE.png');
            imageFilename = '/Users/lopezd/Documents/Babylab/Files/London Eye Tracking Data/Images Eye Tracking/Chair5.png'; 
            PlotFixationsDistance(participantList(participantNumber).fixationX',participantList(participantNumber).fixationY',param.radius,imageFilename,plotFilename)
            %}
    else
    	result=[];
        cellToWrite{participantNumber+1,1} = listParticipants(participantNumber).name;
        cellToWrite{participantNumber+1,2} = 'NaN';
        cellToWrite{participantNumber+1,3} = 'NaN';
        cellToWrite{participantNumber+1,4} = 'NaN';
        cellToWrite{participantNumber+1,5} = 'NaN';
        cellToWrite{participantNumber+1,6} = 'NaN';
        cellToWrite{participantNumber+1,7} = 'NaN';
        cellToWrite{participantNumber+1,8} = 'NaN';
        cellToWrite{participantNumber+1,9} = 'NaN';
        cellToWrite{participantNumber+1,10} = 'NaN';
        cellToWrite{participantNumber+1,11} = 'NaN';
        cellToWrite{participantNumber+1,12} = 'NaN';
        cellToWrite{participantNumber+1,13} = 'NaN';
        cellToWrite{participantNumber+1,14} = 'NaN';
        cellToWrite{participantNumber+1,15} = 'NaN';
        cellToWrite{participantNumber+1,16} = 'NaN';
        cellToWrite{participantNumber+1,17} = 'NaN';
	end
    clear fixations;
end

%Save the data
fid = fopen('Chair6_FixationsRQA_LL4.xls','w');
for listNumber = 1:size(cellToWrite,1)
    if (listNumber == 1)
        for listNumber2 = 1:size(cellToWrite,2)
        	fprintf(fid,'%s\t',cellToWrite{listNumber,listNumber2});
        end
    else
    %fprintf(fid,'%s\t',cellToWrite{listRadius,listNumber,1});
    	for listNumber2 = 1:size(cellToWrite,2)
            if isempty((cellToWrite{listNumber,listNumber2}))
                fprintf(fid,'%s\t','NaN');                    
            elseif isnan(str2num(cellToWrite{listNumber,listNumber2}))
                fprintf(fid,'%s\t',cellToWrite{listNumber,listNumber2});
            else
            	fprintf(fid,'%f\t',str2num(cellToWrite{listNumber,listNumber2}));
            end
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);
