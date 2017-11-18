function runPath = setUpTestFolder(testrunPath, testname)
    mkdir(testrunPath, testname)
    runPath = join([testrunPath, '/', testname]);
end

