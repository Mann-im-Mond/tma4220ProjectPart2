addpath('../utils/')
addpath('../helperFunctions/')

plotter=Plotter(u,mesh,timeInterval);
%Try the following:
plotter.animateSlicePlot2D();
close();
plotter.animateSlicePlot2D(0.01);
close();
plotter.animateSlicePlot2D(0.001);
close();
R=[0,1,0;0,0,1;1,0,0];
v_0 = [0,0,0.1]';
plotter.animateSlicePlot2D(0.01,R,v_0);
close();
plotter.animateSlicePlot3D();
close();
plotter.animateSlicePlot3D(0.01);
close();
plotter.animateSlicePlot3D(0.01,R,v_0);
R=[0,1,0;0,0,1;1,0,0];
v_0 = [0,0,0]';
plotter.animateSlicePlot3D(0.01,R,v_0);
close();
R=[1,-1,0;0,1,0;0,0,1];
v_0 = [0,0,0.075]';
plotter.animateSlicePlot3D(0.01,R,v_0);
close();
plotter.animateScatterPlot();
close();
plotter.shrinkingPlot(452);
close();
%And this is how you can write it into a video file:
for k=1:length(plotter.movies)
    video=VideoWriter(['scliceAnimation',num2str(k),'.avi'],'Uncompressed AVI');
    video.FrameRate=5;
    open(video);
    writeVideo(video,plotter.movies(k).movie)
    close(video);
end
    