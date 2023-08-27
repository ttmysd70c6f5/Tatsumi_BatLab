function RestPosClusteringAll()

original_dir = pwd;
RecFile = dir('*\Ext_Behavior*');
for i = 1:length(RecFile)
    cd(RecFile(i).folder)
    RestPosClustering()
end
cd(original_dir)