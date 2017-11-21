folder  = '../testruns';
list    = dir(folder);
nFolder = length(list);
for runFolder = list'
    runName = runFolder.name;
    if(runName == '.')
        continue
    end
    fullPath = [folder '/' runName];
    stepSize = .1;
    if(contains(runName,'smallSteps'))
        stepSize = .01;
    end
    for dim = 2
        aviName = [runName '_' num2str(dim) 'D.avi'];
        aviPath = [fullPath '/' aviName];
        plotFromFile(fullPath,7200,1,dim,aviPath);
    end
end