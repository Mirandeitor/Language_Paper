function [result,thresholdedDiagonals,thresholdedVerticals]=Cra(x,y,param)
% result=Cra(x,y,param)
%
% parameter:
%
% x                    x 
% y                    y 
% param.delay          e.g. 1
% param.embed          e.g. 3
% param.rescale        0: no rescaling, 2: scaling to mean
% param.metric         'min','max','euclidian','adjacency'
% param.linelength     threshold linelength
% param.radius         e.g. 0.1, ignored with adjacency metric
% param.adjacency      adjacency matrix a(i,j)==1 iff i and j neighbors
%
% result:
%
% result.n             n (length of x,y)
% result.nseq          nseq (number of sequences)
%                      [nseq=n iff delay=1 and embed=1)
% result.xseq          x sequences
% result.yseq          y sequences
% result.recmat        recurrence matrix
% result.nrec          number of recurrences
% result.rec           percentage of recurrences
% result.rr            recurrence rate
% result.det           determinism
% result.meanline      meanline
% result.maxline       maxline
% result.ent           entropy
% result.relent        relative entropy
% result.lam           laminarity
% result.tt            TT
% result.corm          center of recurrence mass
%
% Author               Walter F. Bischof
% Date                 September 2011 - June 2012

% V1.1 The way the distance matrix is calculated using the euclidean
% algorithm has been optimised by David Lopez Perez 21.07.2017

% V 1.2 Bug Fix. There was a problem when the embedded dimension was larger
% than 1 by David Lopez Perez 11.09.2019



    % Default return values
    
    result.n=NaN;
    result.nseq=NaN;
    result.xseq=NaN;
    result.distance=[];
    result.recmat=[];
    result.nrec=NaN;
    result.rec=NaN;
    result.rr=[];
    result.det=NaN;
    result.meanline=NaN;
    result.maxline=NaN;
    result.ent=NaN;
    result.relent=NaN;
    result.lam=NaN;
    result.tt=NaN;
    result.corm=NaN;

    if strcmp(param.metric,'adjacency')
        param.radius=0.5;
    end

    % Ensure that x and y are column vectors (or columns matrices)

    if size(x,1)<size(x,2)
        x=x';
        y=y';
    end
    n=size(x,1);
    m=size(x,2);

    % Are they the same length
    if size(x,1) ~= size(y,1)
        %If not crop the longest one
        if size(x,1) > size(y,1)
            y(size(x,1):end,1) = [];
        else
            x(size(y,1):end,1) = [];
        end
    end
    
    % Compute number of subsequences

    spanSeq=(param.embed-1)*param.delay+1;
    seq=1:param.delay:spanSeq;
    nSeq=n-spanSeq+1;
    nSeq2=nSeq*nSeq;
    lenSeq=length(seq);
    result.n=n;
    result.nseq=nSeq;
    
    if nSeq<=1
        return
    end
    
    % Compute all subsequences
    if ~(lenSeq==1)
        xSeq=zeros(nSeq,lenSeq,m);
        ySeq=zeros(nSeq,lenSeq,m);
    else
        xSeq=zeros(nSeq,m);
        ySeq=zeros(nSeq,m);
    end
    for i=1:nSeq
        xSeq(i,:)=x(i-1+seq,:);
        ySeq(i,:)=y(i-1+seq,:);
    end
    %result.xseq=xSeq;
    %result.yseq=ySeq;
    
    % Compute distance matrix. This uses most of the computation time.

    distanceMatrix = ComputeDistanceMatrix(nSeq,xSeq,ySeq,param.metric,param.adjacency);
    result.distance = distanceMatrix;
    % Rescale (to mean)

    if param.rescale==2
        dmMean=nansum(nansum(distanceMatrix))/nSeq2;
        distanceMatrix=distanceMatrix/dmMean;
        result.distance = distanceMatrix;
    end

    % Recurrence matrix

    recurrenceMatrix=distanceMatrix<=param.radius;
    result.recmat=recurrenceMatrix;
    
    % Diagonal not needed for most measures
    partialRecurrenceMatrix=triu(recurrenceMatrix,1);
    
    % Diagonals of recurrence matrix, excluding main diagonal

    [recurrenceDiagonals,diagonals]=spdiags(partialRecurrenceMatrix);

    % Recurrence measures

    % Global measures: %REC
    nRec=sum(sum(recurrenceMatrix));
    result.nrec=nRec;
    result.rec=100.0*nRec/nSeq2;

    % Recurrence rate
    result.rr=zeros(2,2*nSeq-1);
    result.rr(1,:)=[-nSeq+1:nSeq-1];
    result.rr(2,diagonals+nSeq)=sum(recurrenceDiagonals,1);
    % result.rr(2,:)=result.rr(2,:)/result.rr(1,:);

    % Diagonal line measures: %DET,meanLine,maxLine,TND,ENT
    
    nDiagonals=length(diagonals);
    thresholdedDiagonals=[];
    for i=1:nDiagonals
        d=[0; recurrenceDiagonals(:,i); 0];
        dchange=diff(d);
        dlength=find(dchange==-1)-find(dchange==1);
        dlength=dlength(dlength>=param.linelength);
        thresholdedDiagonals=[thresholdedDiagonals; dlength];
    end

    %Return Diagonals
        
    if ~isempty(thresholdedDiagonals)
        result.det=100*sum(thresholdedDiagonals)/nRec;
        result.meanline=mean(thresholdedDiagonals);
        result.maxline=max(thresholdedDiagonals);
        [result.ent,result.relent]=Entropy(thresholdedDiagonals,param.linelength,result.maxline);
        result.corm=100*sum(sum(recurrenceDiagonals)'.*diagonals)/((nSeq-1)*nRec);
    end
    
    % Horizontal+vertical line measures: %LAM, TT
    % #(verticals + horizontals) in upper triangle = 
    % #(verticals) in upper + #(verticals) in lower triangle
    % Set diagonal to zero to avoid counting over the diagonal
    
    thresholdedVerticals=[];
    paddedRM=[zeros(nSeq,1) recurrenceMatrix-eye(nSeq) zeros(nSeq,1)];
    % In CRA, diagonal elements can be zero 
    paddedRM=max(paddedRM,0);
    for i=1:nSeq
        v=paddedRM(i,:);
        vchange=diff(v);
        vlength=find(vchange==-1)-find(vchange==1);
        vlength=vlength(vlength>=param.linelength);
        thresholdedVerticals=[thresholdedVerticals vlength];
    end
    
    %Return Vertical Lines
    
    if ~isempty(thresholdedVerticals)
        result.lam=100*sum(thresholdedVerticals)/sum(sum(recurrenceMatrix));
        result.tt=mean(thresholdedVerticals);
    end
end

function [entropy,relativeEntropy]=Entropy(a,minLine,maxLine)
    [p,~]=hist(a,unique(a));
    p=p/sum(p);
    entropy=-sum(p.*log2(p));
    if maxLine==minLine
        relativeEntropy=NaN;
    else
        relativeEntropy=entropy/log2(maxLine-minLine+1);
    end
end

function distanceMatrix=ComputeDistanceMatrix(nSeq,xSeq,ySeq,metric,adjacency)
% Compute the distance matrix between xSeq and ySeq.
% This routine uses half the computation time and needs optimization,
% probably along the lines of distance.m. xSeq and ySeq are of the form
% xSeq(n,M,2) n = number of sequences, M = embedding dimension, last for
% (x,y)

    distanceMatrix=zeros(nSeq,nSeq);
    if strcmp(metric,'euclidian')
        %{
        for i=1:nSeq
            for j=1:nSeq
                distanceMatrix(i,j)=sqrt(sum((xSeq(i,:)-ySeq(j,:)).^2,2));
            end
        end
        %}
        distanceMatrix = pdist2(reshape(xSeq,[nSeq size(xSeq,2)]) ,reshape(ySeq,[nSeq size(ySeq,2)]));
    elseif strcmp(metric,'max')
        for i=1:nSeq
            for j=1:nSeq
                distanceMatrix(i,j)=max(abs(xSeq(i,:)-ySeq(j,:)),[],2);
            end
        end
    elseif strcmp(metric,'min')
        for i=1:nSeq
            for j=1:nSeq
                distanceMatrix(i,j)=min(abs(xSeq(i,:)-ySeq(j,:)),[],2);
            end
        end
    elseif strcmp(metric,'adjacency')
        for i=1:nSeq
            for j=1:nSeq
                % We must deal with outliers, which are not adjacent to anything
                xx=xSeq(i,:);
                yy=ySeq(j,:);
                d=ones(1,length(xx));
                okValues=find(xx>=0 & yy>=0);
                d(okValues)=1-adjacency(xx(okValues),yy(okValues));
                distanceMatrix(i,j)=d;
            end
        end
    else
        disp(['Unknown distance metric:' metric]);
    end
end
