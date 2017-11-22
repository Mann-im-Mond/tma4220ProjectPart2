classdef Plotter < handle
% Plotter(u,mesh,timeInterval) generates an object, that can be used to
% generate and store several plots of the function u on a mesh and animated
% plots over the timeInterval.
%
% u:        NxM Matrix of a descrete function, where N is equal to the
%           number of points in mesh and M is equal to the stored time
%           steps from timeInterval
% mesh:     FullMesh Object, that includes the space descretization.
% timeInterval:
%           TimeInterval Object, that includes the time descretization.
% dim:      Dimension of the mesh to be plottet.
% movies:   Struct containg created movies, that can be called with
%           movie(movies(i).movie)
% u_min:    Value that is smaller or at most equal to the minum over all u
% u_min:    Value that is greater or at least equal to the maximum over all u
% 
% See also FullMesh,TimeInterval

    properties
        u
        mesh
        timeInterval
        dim
        movies
        u_min
        u_max
    end
    
    methods (Access=public,Static)
        function rgb=colorConverterStatic(u,u_min,u_max)
        %colorConverterStatic(u,u_min,u_max) converts the given vector u
        %into a color vector, where each entry of u will be converted into
        %an rgb triplet where u_min will be converted to blue, u_max will
        %be conderted to red and the middle (u_max-u-min)/2 will be
        %converted to yellow inbetween the colors will be linear
        %interpolated. So entries of u have to be between u_min and u_max.
            %Shrink u into the interval [0,1]
            to01Convert=@(v) (v-u_min)/(u_max-u_min);
            %Create RGB triplet.
            cC=@(v) [min(2*to01Convert(v),1),min(2*to01Convert(v),2-2*to01Convert(v)),max(1-2*to01Convert(v),0)];
            rgb=cC(u);
        end
    end
    
    methods (Access=private)
        function appendMovie(obj,movie)
        %appendMovie appends a given movie to the movies struct.
            [~,tmp]=size(obj.movies);
            obj.movies(tmp+1).movie=movie;
        end
    end
    
    methods
        function rgb=colorConverter(obj,K)
        %colorConverter(K) calls colorConverterStatic for time step K and
        %u_min and u_max from the object and returns a vector with rgb
        %triplets.
        %
        %See also colorConverterStatic
            rgb=Plotter.colorConverterStatic(obj.u(:,K),obj.u_min,obj.u_max);
        end
        
        function obj=Plotter(u,mesh,timeInterval)
            %Standard constructor
            [N,M]=size(u);
            if not(N==length(mesh.points(:,1)))
                error('The number of rows must be equal to the number of points in mesh');
            end
            if (M<1)
                error('u must not be empty');
            elseif (M==1)
                timeInterval.t_max=timeInterval.t_0;
                timeInterval.setDescreteInterval();
            elseif (M<timeInterval.getNumberOfSteps()/timeInterval.n_to_plot+1)
                timeInterval.t_max=(M-1)*timeInterval.h;
                timeInterval.setDescreteInterval();
            elseif (M>timeInterval.getNumberOfSteps()/timeInterval.n_to_plot+1)
                error('To many timesteps');
            end
            obj.u=u;
            obj.mesh=mesh;
            obj.timeInterval=timeInterval;
            [~,obj.dim]=size(mesh.points);
            obj.u_min=min(min(u));
            obj.u_max=max(max(u));
            obj.movies=[];
        end
        
        function addLegend(obj,K)
        %addLegend(K) adds a legend for the k-th step to the current figure.
            legend(['t=',num2str(obj.timeInterval.descreteInterval((K-1)*obj.timeInterval.n_to_plot+1)/60) 'min' newline...
                'u_{average}=' num2str(mean(obj.u(:,K))-273) '°C' newline...
                'u_{min}=' num2str(min(obj.u(:,K))-273) '°C' newline...
                'u_{max}=' num2str(max(obj.u(:,K))-273) '°C'],'Location','northoutside');
        end
        
        function addLabel(obj)
        %addLabel() adds axis lables to the present axis.
            if(obj.dim==2)
                xlabel([num2str(min(obj.mesh.points(:,1))),'<= x <= ',num2str(max(obj.mesh.points(:,1)))]); % x-axis label
                ylabel([num2str(min(obj.mesh.points(:,2))),'<= y <= ',num2str(max(obj.mesh.points(:,2)))]); % y-axis label
            elseif(obj.dim==3)
                xlabel([num2str(min(obj.mesh.points(:,1))),'<= x <= ',num2str(max(obj.mesh.points(:,1)))]); % x-axis label
                ylabel([num2str(min(obj.mesh.points(:,2))),'<= y <= ',num2str(max(obj.mesh.points(:,2)))]); % y-axis label
                zlabel([num2str(min(obj.mesh.points(:,3))),'<= z <= ',num2str(max(obj.mesh.points(:,3)))]); % z-axis label
            else
                error('Can only plot 2 or 3 dimensions.');
            end
        end
        
        function fig1=initAnimation(obj,name)
        %initAnimation(name) creates a figure and gives it a title as well
        %as a Label for the axis. It also makes sure, that for the 3D slice
        %plot the plot will have three dimensions.
            fig1=figure('Name',[name,' plot animation']);
            title([name,' plot animation']);
            if(strcmp(name,'3D slice'))
                view(3);
            end
            obj.addLabel();
        end
        
        function movie=animatePlot(obj,name,varargin)
        %animatePlot(name,varargin) will call the function specified in
        %name with values from varargin and return a movie.
            switch name
                case 'Scatter'
                    plotterFunction=@(K) obj.scatterPlotSingleTimeStep(K);
                case 'Shrinking'
                    threshold=varargin{1};
                    plotterFunction=@(K) obj.shrinkingPlotSingleTimeStep(threshold,K);
                case 'Slice2D'
                    delTri=varargin{1};
                    plainPoints=varargin{2};
                    plotterFunction=@(K) obj.slicePlotSingleTimeStep2D(delTri,plainPoints,K);
                    name='2D slice';
                case 'Slice3D'
                    delTri=varargin{1};
                    plainPoints=varargin{2};
                    plotterFunction=@(K) obj.slicePlotSingleTimeStep3D(delTri,plainPoints,K);
                    name='3D slice';
                otherwise
                    error([name,' not recognized']);
            end
            %Initialize the animation
            fig1=obj.initAnimation(name);
            %Plot the first step
            plotterFunction(1);
            %Fix axis
            axis tight manual;
            ax = gca;
            ax.NextPlot = 'replaceChildren';
            drawnow;
            %Get first frame
            movie(1) = getframe(fig1);
            [~,loops]=size(obj.u);
            movie(loops) = struct('cdata',[],'colormap',[]);
            %Loop over all time steps and store frames
            for j = 2:loops
                plotterFunction(j);
                drawnow;
                movie(j) = getframe(fig1);
            end
            %Append the movie to movies
            obj.appendMovie(movie);
        end
        
        %------------------Functions for scatter plot---------------------
        
        function handler=scatterPlotSingleTimeStep(obj,K)
        %scatterPlotSingleTimeStep(K) calls scatterPlotSingleTimeStep2D if
        %dim=2 and scatterPlotSingleTimeStep3D if dim=3.
            if(obj.dim==2)
                handler=obj.scatterPlotSingleTimeStep2D(K);
            elseif(obj.dim==3)
                handler=obj.scatterPlotSingleTimeStep3D(K);
            else
                error('Only dimension 2 and 3 can be plotted.');
            end
        end
        
        function handler=scatterPlotSingleTimeStep2D(obj,K)
        %scatterPlotSingleTimeStep(K) plots a heat scatter plot in two
        %dimensions
            if(not(obj.dim==2))
                error('The dimension of the mesh should be 2');
            end
            if (K<1 || K>length(obj.u(1,:)))
                error(['K has to be between 1 and ',length(obj.u(1,:))]);
            end
            handler=scatter(obj.mesh.points(:,1),obj.mesh.points(:,2),10.*ones(length(obj.mesh.points(:,1)),1),obj.colorConverter(K),'filled');
            obj.addLegend(K);
        end
        
        function handler=scatterPlotSingleTimeStep3D(obj,K)
        %scatterPlotSingleTimeStep(K) plots a heat scatter plot in three
        %dimensions
            if(not(obj.dim==3))
                error('The dimension of the mesh should be 3');
            end
            if (K<1 || K>length(obj.u(1,:)))
                error(['K has to be between 1 and ',length(obj.u(1,:))]);
            end
            handler=scatter3(obj.mesh.points(:,1),obj.mesh.points(:,2),obj.mesh.points(:,3),10.*ones(length(obj.mesh.points(:,1)),1),obj.colorConverter(K),'filled');
            obj.addLegend(K);
        end
                    
        function movie=animateScatterPlot(obj)
        %animateScatterPlot() will create an animated scatter heat plot.
            movie=obj.animatePlot('Scatter');
        end
        
        %-----------------Functions for shrinking plot--------------------
        
        function handler=shrinkingPlotSingleTimeStep(obj,threshold,K)
        %shrinkingPlotSingleTimeStep(threshold,K)Plots the boundary of the
        %area/volume in which u<=threshold in timestep k. This is realy
        %slow and does not look as good as hoped
            %Find the position of the values of u where u<=threshold
            position=find(obj.u(:,K)<=threshold);
            %If all points fullfill u<=threshold plot empty plot
            if(isempty(position))
                handler=obj.boundaryPlot([]);
                obj.addLegend(K);
                return;
            end
            %Create vector, that projects to the remaining points
            projection=zeros(length(obj.u(:,K)),1);
            counter=1;
            for k=1:length(obj.u(:,K))
                if(counter<=length(position(:,1)) && position(counter)==k)
                    projection(k)=counter;
                    counter=counter+1;
                end
            end
            %Set the remaining points
            points=obj.mesh.points(position,:);
            %Set the remaining triangles 
            [i,j]=size(obj.mesh.triangulation);
            triangles=zeros(i,j);
            counter=0;
            for k=1:i
                if(all(ismember(obj.mesh.points(obj.mesh.triangulation(k,:),:), points,'rows')))
                    counter=counter+1;
                    triangles(counter,:)=projection(obj.mesh.triangulation(k,:));
                end
            end
            triangles=triangles(1:counter,:);
            %remove loose points
            projection=zeros(length(points(:,1)),1);
            counter=0;
            for k=1:length(points(:,1))
                if(ismember(k,triangles))
                    projection(k)=k-counter;
                    points(k-counter)=points(k);
                else
                    counter=counter+1;
                end
            end
            points=points(1:k-counter,:);
            for k=1:length(triangles(:,1))
                triangles(k,:)=projection(triangles(k,:));
            end
            %Create triangulation
            if(length(points(:,1))<=4)
                handler=obj.boundaryPlot([]);
                obj.addLegend(K);
                return
            end
            tri=triangulation(triangles,points);
            %Find boundry of the triangulation and plot it
            handler=obj.boundaryPlot(tri.freeBoundary());
            obj.addLegend(K);
        end
        
        function handler=boundaryPlot(obj,boundary)
        %boundaryPlot(boundary) calls boundaryPlot2D(boundary) if dim==2
        %and boundaryPlot3D(boundary) if dim==3.
            if(obj.dim==2)
                handler=obj.boundaryPlot2D(boundary);
            elseif(obj.dim==3)
                handler=obj.boundaryPlot3D(boundary);
            else
                error('Only dimension 2 and 3 can be plotted.');
            end
        end
        
        function handler=boundaryPlot2D(obj,boundary)
        %boundaryPlot2D(boundary) plots the given boundary within the
        %2D-mesh in obj. The boundary does not need to be the boundary of
        %the mesh.
            %Find path through the points that represent the boundary
            if(isempty(boundary))
                %Plot nothing
                handler=scatter3(obj.mesh.points(1,1),obj.mesh.points(1,2),obj.mesh.points(1,3),'MarkerEdgeColor','white');
                legend(['t=',num2str(obj.timeInterval.descreteInterval((K-1)*obj.timeInterval.n_to_plot+1))],'Location','northoutside');
                return;
            end
            faces=zeros(1,length(boundary(:,1)));
            faces(1)=boundary(1,1);
            positionX=1;
            positionY=1;
            for k=2:length(boundary(:,1))
                faces(k)=boundary(positionY,3-positionX);
                next=find(boundary(:,positionX)==faces(k));
                if(isempty(next))
                    positionX=3-positionX;
                    next=find(boundary(:,positionX)==faces(k));
                    positionY=next(next~=positionY);
                elseif(length(next(:,1))==1)
                    positionY=next;
                else
                    error('The boundary is not a circle.');
                end
            end
            if(faces(1)~=boundary(positionY,3-positionX))
                error('The boundary is not a circle.');
            end
            %plot the boundary
            handler=patch('Faces',faces,'Vertices',obj.mesh.points,'EdgeColor','blue','FaceColor','none');
        end
        
        function handler=boundaryPlot3D(obj,boundary)
        %boundaryPlot3D(boundary) plots the given boundary within the
        %3D-mesh in obj. The boundary does not need to be the boundary of
        %the mesh.
            if(isempty(boundary))
                %Plot nothing
                handler=scatter3(obj.mesh.points(1,1),obj.mesh.points(1,2),obj.mesh.points(1,3),'MarkerEdgeColor','white');
                return;
            else
                l=length(boundary(:,1));
            end
            %Initialize
            X=zeros(3,l);
            Y=zeros(3,l);
            Z=zeros(3,l);
            C=zeros(3,l);
            %Set the color to blue everywhere
            C(:,3)=1;
            for k=1:l
                X(:,k)=obj.mesh.points(boundary(k,:),1);
                Y(:,k)=obj.mesh.points(boundary(k,:),2);
                Z(:,k)=obj.mesh.points(boundary(k,:),3);
            end
            %plot the boundary
            handler=fill3(X,Y,Z,C);
        end
        
        function movie=shrinkingPlot(obj,threshold)
        %shrinkingPlot(threshold) plots an animation over the timesteps in
        %u where u<=threshold. The animation will plot the bounary of the
        %points where u<=thresholds holds.
            movie=obj.animatePlot('Shrinking',threshold);
        end
        
        %--------------------Functions for slice plots--------------------
        
        function movie=animateSlicePlot2D(obj,varargin)
        %animateSlicePlot2D(vargin) Plots a heat animation on a slice
        %plane over the timesteps in u. The slice plane is by default the
        %standard plane where x=0, but can also be given in varargin. The
        %density of points is by default given by 0.005
        %Options for varargin:
        %varargin=(density)
        %       where density is the fineness of the plot
        %varargin=(R,v_o)
        %       where the plain is given by R, the rotation matrix of
        %       compared to the standard plane and v_0 is the moving
        %varargin=(density,R,v_o)
        %       where density,R and v_o have the same meaning as above
            if(nargin==1)
                density=0.005;
                plain=Plain(eye(3),zeros(3,1));
            elseif(nargin==2)
                density=varargin{1};
                plain=Plain(eye(3),zeros(3,1));
            elseif(nargin==3)
                density=0.005;
                plain=Plain(varargin{1},varargin{2});
            elseif(nargin==4)
                density=varargin{1};
                plain=Plain(varargin{2},varargin{3});
            else
                error('animateSlicePlot needs 0,1,2 or 3 input values (next to the plotter object)');
            end
            delTri=getPlainTriangulation(obj,plain,density);
            plainPoints=zeros(length(delTri.Points(:,1)),3);
            for k=1:length(delTri.Points(:,1))
                plainPoints(k,:)=plain.plain(delTri.Points(k,1),delTri.Points(k,2))';
            end
            
            movie=obj.animatePlot('Slice2D',delTri,plainPoints);
        end
        
        function handler=slicePlotSingleTimeStep2D(obj,delTri,plainPoints,K)
            interpol = scatteredInterpolant(obj.mesh.points,obj.u(:,K));
            u_test=interpol(plainPoints(:,1),plainPoints(:,2),plainPoints(:,3));
            C=Plotter.colorConverterStatic(u_test,obj.u_min,obj.u_max);
            handler=patch('Faces',delTri.ConnectivityList,'Vertices',delTri.Points,'FaceColor','flat','FaceVertexCData',C,'EdgeColor','none');
            obj.addLegend(K);
        end
        
        function movie=animateSlicePlot3D(obj,varargin)
        %Plots an animation over the timesteps in u on the slice where x=0.
            if(nargin==1)
                density=0.005;
                plain=Plain(eye(3),zeros(3,1));
            elseif(nargin==2)
                density=varargin{1};
                plain=Plain(eye(3),zeros(3,1));
            elseif(nargin==3)
                density=0.005;
                plain=Plain(varargin{1},varargin{2});
            elseif(nargin==4)
                density=varargin{1};
                plain=Plain(varargin{2},varargin{3});
            else
                error('animateSlicePlot needs 0,1,2 or 3 input values (next to the plotter object)');
            end
            delTri=getPlainTriangulation(obj,plain,density);
            plainPoints=zeros(length(delTri.Points(:,1)),3);
            for k=1:length(delTri.Points(:,1))
                plainPoints(k,:)=plain.plain(delTri.Points(k,1),delTri.Points(k,2))';
            end
            
            movie=obj.animatePlot('Slice3D',delTri,plainPoints);
        end
        
        function handler=slicePlotSingleTimeStep3D(obj,delTri,plainPoints,K)
            interpol = scatteredInterpolant(obj.mesh.points,obj.u(:,K));
            u_test=interpol(plainPoints(:,1),plainPoints(:,2),plainPoints(:,3));
            C=Plotter.colorConverterStatic(u_test,obj.u_min,obj.u_max);
            handler=patch('Faces',delTri.ConnectivityList,'Vertices',plainPoints,'FaceColor','flat','FaceVertexCData',C,'EdgeColor','none');
            obj.addLegend(K);
        end
        
        function delTri=getPlainTriangulation(obj,plain,density)
            %Rotate
            points2D=plain.plain_inverse(obj.mesh.points')';
            %Project to plain
            points2D=points2D(:,2:3);
            %Get the boundary
            boundary2D=boundary(points2D);
            boundarypolygon=points2D(boundary2D,:);
            %Create square mesh
            X=(min(points2D(:,1)):density:max(points2D(:,1)))';
            Y=(min(points2D(:,2)):density:max(points2D(:,2)))';
            square_mesh=zeros(length(X)*length(Y),2);
            for k=1:length(X)
                square_mesh((k-1)*length(Y)+1:k*length(Y),:)=[X(k)*ones(length(Y),1),Y];
            end
            %create mesh within boundary
            delTri=delaunayTriangulation(square_mesh(inpolygon(square_mesh(:,1),square_mesh(:,2),boundarypolygon(:,1),boundarypolygon(:,2)),:));
        end
    end
end