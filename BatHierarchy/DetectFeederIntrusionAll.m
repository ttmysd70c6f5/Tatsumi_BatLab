function DetectFeederIntrusionAll()

original_dir = pwd;
RecFile = dir('*\Ext_Behavior*');
for i = 1:length(RecFile)
    cd(RecFile(i).folder)
    DetectFeederIntrusion
    pause(1)
end
cd(original_dir)