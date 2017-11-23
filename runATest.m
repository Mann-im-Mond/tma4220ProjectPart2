addpath('utils/')
addpath('helperFunctions/')

neumannIdentifier = @(x) false;
mesh = getCakeMesh('grid', neumannIdentifier, .15);
alpha = @(isRod) alphaFkt(isRod);

t_max = 7200;
timeStep = 6;
plot_num = 60/timeStep;
timeInterval = TimeInterval(0,t_max,timeStep,plot_num);
u0 = @(x) x(1)*0 +293;
gD = @(x,t) x(1) + 473;
gN = @(x,t) x(1)*0;
solver = HeatEquationSolver(mesh,timeInterval,alpha,u0,gD,gN,@(x) false);
u = solver.solve();
plotter=Plotter(u,mesh,timeInterval);
plotter.animateSlicePlot2D();
close();
plotter.animateSlicePlot2D([0,1,0;0,0,1;1,0,0],[0,0,0.1]');
close();
plotter.animateSlicePlot3D([0,1,0;0,0,1;1,0,0],[0,0,0.1]');
close();
plotter.animateScatterPlot();
close();
diary off