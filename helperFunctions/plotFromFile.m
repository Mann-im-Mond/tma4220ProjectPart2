function [] = plotFromFile(folderPath,t_max,timeStep,dim)
    u=dlmread([folderPath '/temperatureSolution.txt']);
    mesh = getCakeMesh('../grid', @(x) false, .15);
    plot_num = 60/timeStep;
    timeInterval = TimeInterval(0,t_max,timeStep,plot_num);
    plotter=Plotter(u,mesh,timeInterval);
    if(dim == 2)
        plotter.animateSlicePlot2D();
    else if(dim == 3)
        plotter.animateScatterPlot();
    end
end

