function TrjCluster = adhocMergeTrajectory(TrjCluster,FigDir,recName,expType)

if all([strcmp(recName,'240207'),strcmp(expType,'social')])
    TrjCluster = mergeTrajectoryCluster(TrjCluster,FigDir,6,8,12); % bat id = 6, merge cluster 8 and 12
    TrjCluster = mergeTrajectoryCluster(TrjCluster,FigDir,7,4,6); % bat id = 6, merge cluster 8 and 12
    TrjCluster = mergeTrajectoryCluster(TrjCluster,FigDir,7,22,16); % bat id = 6, merge cluster 8 and 12
elseif all([strcmp(recName,'240208'),strcmp(expType,'social')])
    TrjCluster = mergeTrajectoryCluster(TrjCluster,FigDir,6,2,23,50); % bat id = 6, merge cluster 8 and 12
    TrjCluster = mergeTrajectoryCluster(TrjCluster,FigDir,6,25,4,51); % bat id = 6, merge cluster 8 and 12
    TrjCluster = mergeTrajectoryCluster(TrjCluster,FigDir,6,1,8,52); % bat id = 6, merge cluster 8 and 12
    TrjCluster = mergeTrajectoryCluster(TrjCluster,FigDir,6,52,22,53); % bat id = 6, merge cluster 8 and 12
    TrjCluster = mergeTrajectoryCluster(TrjCluster,FigDir,7,5,19,50); % bat id = 6, merge cluster 8 and 12
    TrjCluster = mergeTrajectoryCluster(TrjCluster,FigDir,7,16,17,51); % bat id = 6, merge cluster 8 and 12
    TrjCluster = mergeTrajectoryCluster(TrjCluster,FigDir,7,1,3,52); % bat id = 6, merge cluster 8 and 12
    TrjCluster = mergeTrajectoryCluster(TrjCluster,FigDir,7,15,52,53); % bat id = 6, merge cluster 8 and 12
end

end