function [] = runBakeCake(neumannIdentifier,timeStep,t_max,gD,gN,stoppingCondition,name)
    addpath('../utils/')
    addpath('../helperFunctions/')

    dateTimeStamp = getDateTimeStamp();
    outputFolder = setUpTestFolder('../testruns',[dateTimeStamp '_' name]);
    diary([outputFolder '/logfile.log'])
    diary on

    dateTimeStamp
    neumannIdentifier
    mesh = getCakeMesh('../grid', neumannIdentifier, .03)
    alpha = @(isRod) alphaFkt(isRod)
    plot_num = 60/timeStep;
    timeInterval = TimeInterval(0,t_max,timeStep,plot_num)
    u0 = @(x) x(1)*0 +293
    gD
    gN
    solver = HeatEquationSolver(mesh,timeInterval,alpha,u0,gD,gN,stoppingCondition);
    %u=dlmread('../testruns/18Nov2017_18-32-01/temperatureSolution.txt');
    disp('start to solve.')
    u = solver.solve();
    %plotter=Plotter(u,mesh,timeInterval);
    %plotter.animateScatterPlot;
    %plotter.shrinkingPlotSingleStep(plotter.u_max,1);
    %plotter.shrinkingPlot(452); %If the cake is done at 452 this makes sense
    %but this will take a lot of time and might not look good.
    dlmwrite([outputFolder '/temperatureSolution.txt'],u);
    diary off
end

