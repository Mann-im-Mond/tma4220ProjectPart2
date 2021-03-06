function movieOut = plotFromFile(folderPath,t_max,timeStep,dim,outputFile)
    u=dlmread([folderPath '/temperatureSolution.txt']);
    mesh = getCakeMesh('../grid', @(x) false, .15);
    plot_num = 60/timeStep;
    timeInterval = TimeInterval(0,t_max,timeStep,plot_num);
    plotter=Plotter(u,mesh,timeInterval);
    if(dim == 2)
        movieOut = plotter.animateSlicePlot2D(0.001);
    elseif(dim == 2.4)
        movieOut = plotter.animateSlicePlot2D(0.001,[0,1,0;0,0,1;1,0,0],[0,0,0.1]');
    elseif(dim == 2.6)
        movieOut = plotter.animateSlicePlot3D(0.001,[0,1,0;0,0,1;1,0,0],[0,0,0.1]');
    elseif(dim == 3)
        movieOut = plotter.animateScatterPlot();
    end
    close();
    if not(isempty(outputFile))
        disp(['writing movie to ' outputFile])
        v = VideoWriter(outputFile);
        v.FrameRate = 6;
        open(v);
        writeVideo(v,movieOut)
        close(v)
    end
end