addpath('../utils/')
addpath('../helperFunctions/')

dateTimeStamp = getDateTimeStamp();
outputFolder = setUpTestFolder('../testruns',[dateTimeStamp '_alpha']);
diary([outputFolder '/logfile.log'])
diary on

dateTimeStamp
neumannIdentifier = @(x) (x(1)*0)
mesh = getCakeMesh('../grid', neumannIdentifier, .15)
alpha = @(isRod) alphaFkt(isRod)
timeInterval = TimeInterval(0,3600,1,60)
u0 = @(x) x(1)*0 +293
gD = @(x,t) norm(x(:)'*[0;0;1],1)*0 +453 %-180*exp(-3.5*t/600)
gN = @(x,t) norm(x(:)'*[1;0;1],1)*0 +1 +0*t
solver = HeatEquationSolver(mesh,timeInterval,alpha,u0,gD,gN);
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