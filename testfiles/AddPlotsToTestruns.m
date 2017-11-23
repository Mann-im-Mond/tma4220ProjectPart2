folder  = '../testruns';
list    = dir(folder);
nFolder = length(list);
for runFolder = list'
    runName = runFolder.name;
    if(runName == '.')
        continue
    end
    if not(contains(runName, 'pureNormalCake'))
        continue
    end
    fullPath = [folder '/' runName];
    t_max = 7200;
    stepSize = .1;
    if(contains(runName,'smallSteps'))
        stepSize = .01;
    end
    if contains(runName, 'pure')
        t_max = 14400;
    end
    for dim = [2 2.4 2.6 3]
        aviName = [runName '_' num2str(dim) 'D_fine.avi'];
        aviPath = [fullPath '/' aviName];
        plotFromFile(fullPath,t_max,stepSize,dim,aviPath);
    end
end