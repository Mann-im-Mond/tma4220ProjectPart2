addpath('../utils/')
addpath('../helperFunctions/')

dateTimeStamp = getDateTimeStamp();
outputFolder = setUpTestFolder('../testruns',dateTimeStamp);
diary(join([outputFolder '/logfile.log']))
diary on

dateTimeStamp
neumannIdentifier = @(x) (x(1)*0)
mesh = getCakeMesh('../grid', neumannIdentifier, .15)
alpha = @(x) x(1)*0 + 1e-7
timeInterval = TimeInterval(0,1200,.1,600)
u0 = @(x) x(1)*0 +293
gD = @(x,t) norm(x(:)'*[0;0;1],1)*0 +453 +0*t^2
gN = @(x,t) norm(x(:)'*[1;0;1],1)*0 +1 +0*t
solver = HeatEquationSolver(mesh,timeInterval,alpha,u0,gD,gN);
u = solver.solve();
dlmwrite(join([outputFolder '/temperatureSolution.txt']),u);
diary off