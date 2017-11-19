function runPath = setUpTestFolder(testrunPath, testname)
    mkdir(testrunPath, testname)
    runPath = [testrunPath, '/', testname];
end

