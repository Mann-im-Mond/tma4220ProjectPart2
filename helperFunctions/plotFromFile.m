function movieOut = plotFromFile(folderPath,t_max,timeStep,dim,outputFile)
    u=dlmread([folderPath '/temperatureSolution.txt']);
    mesh = getCakeMesh('../grid', @(x) false, .15);
    plot_num = 60/timeStep;
    timeInterval = TimeInterval(0,t_max,timeStep,plot_num);
    plotter=Plotter(u,mesh,timeInterval);
    if(dim == 2)
        movieOut = plotter.animateSlicePlot2D();
    elseif(dim == 3)
        movieOut = plotter.animateScatterPlot();
    end
    if not(isempty(outputFile))
        disp(['writing movie to ' outputFile])
        v = VideoWriter(outputFile);
        v.FrameRate = 6;
        open(v);
        writeVideo(v,movieOut)
        close(v)
    end
end

