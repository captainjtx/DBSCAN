function demo_stepwise()
%twospirals
data=twospirals(200,360,50,1.5,15);
% data=corners(500);
s=zeros(size(data,1),1);

for i=1:size(data,1)
    hold on
    s(i)=scatter(data(i,1),data(i,2),'filled','markerfacecolor',[0.8,0.8,0.8]);
end
drawnow

%%
%using euclidean distance
distmat=zeros(size(data,1),size(data,1));

for i=1:size(data,1)
    for j=i:size(data,1)
        distmat(i,j)=sqrt((data(i,1:2)-data(j,1:2))*(data(i,1:2)-data(j,1:2))');
    end
end

for i=1:size(data,1)
    for j=i:size(data,1)
        distmat(j,i)=distmat(i,j);
    end
end

% k_dist=zeros(size(data,1),1);
% figure
% for k=3:5
%     for i=1:size(data,1)
%         tmp=sort(distmat(:,i),'ascend');
%         k_dist(i)=tmp(k);
%     end
%     hold on
%     plot(1:size(data,1),k_dist);
% end
%%
Eps=0.5;
MinPts=4;

DBSCAN_STEPWISE(s,distmat,Eps,MinPts);

end

function Clust = DBSCAN_STEPWISE(s,DistMat,Eps,MinPts)

%A step-wise illustration of DBSCAN on 2D data

%A simple DBSCAN implementation of the original paper:
%"A Density-Based Algorithm for Discovering Clusters in Large Spatial
%Databases with Noise" -- Martin Ester et.al.
%Since no spatial access method is implemented, the run time complexity
%will be N^2 rather than N*logN
%**************************************************************************
%Input: DistMat, Eps, MinPts
%DistMat: A N*N distance matrix, the (i,j) element contains the distance
%from point-i to point-j.

%Eps:     A scalar value for Epsilon-neighborhood threshold.

%MinPts:  A scalar value for minimum points in Eps-neighborhood that holds
%the core-point condition.
%**************************************************************************
%Output: Clust
%Clust:  A N*1 vector describes the cluster membership for each point. 0 is
%reserved for NOISE.
%**************************************************************************
%Written by Tianxiao Jiang, jtxinnocence@gmail.com
%Nov-4-2015
%**************************************************************************

%Initialize Cluster membership as -1, which means UNCLASSIFIED
Clust=zeros(size(DistMat,1),1)-1;
ClusterId=1;
ClusterColor=rand(1,3);

%randomly choose the visiting order
VisitSequence=randperm(length(Clust));

for i=1:length(Clust)
    % For each point, check if it is not visited yet (unclassified)
    pt=VisitSequence(i);
    if Clust(pt)==-1
        %Iteratively expand the cluster through density-reachability
        [Clust,isnoise]=ExpandCluster(s,DistMat,pt,ClusterId,Eps,MinPts,Clust,ClusterColor);
        if ~isnoise
            ClusterId=ClusterId+1;
            ClusterColor=rand(1,3);
        end
    end
end
end


function [Clust,isnoise]=ExpandCluster(s,DistMat,pt,ClusterId,Eps,MinPts,Clust,ClusterColor)

%region query
seeds=find(DistMat(:,pt)<=Eps);

if length(seeds)<MinPts
    Clust(pt)=0; % 0 reserved for noise
    set(s(pt),'Marker','*');
    pause(0.01)
    isnoise=true;
    return
else
    Clust(seeds)=ClusterId;
    set(s(seeds),'MarkerFaceColor',ClusterColor);
    pause(0.01)
    %delete the core point
    seeds=setxor(seeds,pt);
    while ~isempty(seeds)
        currentP=seeds(1);
        %region query
        result=find(DistMat(:,currentP)<=Eps);
        if length(result)>=MinPts
            for i=1:length(result)
                resultP=result(i);
                if Clust(resultP)==-1||Clust(resultP)==0 % unclassified or noise
                    set(s(resultP),'MarkerFaceColor',ClusterColor);
                    pause(0.01)
                    if Clust(resultP)==-1 %unclassified
                        seeds=[seeds(:);resultP];
                    end
                    Clust(resultP)=ClusterId;
                end
                
            end
        end
        seeds=setxor(seeds,currentP);
    end
    isnoise=false;
    
    return
end
end

