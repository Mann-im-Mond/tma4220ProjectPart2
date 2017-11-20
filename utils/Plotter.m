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
% movies:   Vector with struct containg created movies, that can be called
%           with movie(movies(i))
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
    
    properties (Access=private)
        pdeModel
    end
    
    methods (Access=private,Static)
        function rgb=colorConverterStatic(u,min,max)
            to01Convert=@(v) (v-min)/(max-min);
            cC=@(v) [min(2*to01Convert(v),1),min(2*to01Convert(v),2-2*to01Convert(v)),max(1-2*to01Convert(v),0)];
            rgb=cC(u);
        end
        
        function handler=slicePlotStatic(u,tri,plain_points)
            l=length(tri);
            Y=zeros(3,l);
            Z=zeros(3,l);
            %Set the color to blue everywhere
            for k=1:l
                Y(:,k)=plain_points(tri(1,:),1);
                Z(:,k)=plain_points(tri(:,k),2);
            end
            C=Plotter.colorConverterStatic(u,min(u),max(u));
            handler=fill3(zeros(1,3),Y,Z,C,'filled');
        end
    end
    
    methods (Access=private)
        function appendMovie(obj,movie)
            [~,tmp]=size(obj.movies);
            obj.movies(tmp+1).movie=movie;
        end
    end
    
    methods
        function rgb=colorConverter(obj,K)
            rgb=Plotter.colorConverterStatic(obj.u(:,K),obj.u_min,obj.u_max);
        end
        function obj=Plotter(u,mesh,timeInterval)
            [N,M]=size(u);
            if not(N==length(mesh.points(:,1)))
                error('The number of rows must be equal to the number of points in mesh');
            end
            if (M<1)
                error('u must not be empty');
            elseif (M==1)
                timeInterval.t_max=timeInterval.t_0;
                timeInterval.setDescreteInterval();
            elseif (M<timeInterval.getNumberOfSteps()+1)
                timeInterval.t_max=(M-1)*timeInterval.h;
                timeInterval.setDescreteInterval();
            end
            obj.u=u;
            obj.mesh=mesh;
            obj.timeInterval=timeInterval;
            [~,obj.dim]=size(mesh.points);
%             pdeModel=createpde();
%             pdeModel.geometryFromMesh(mesh.points',mesh.triangulation');
            obj.u_min=min(min(u));
            obj.u_max=max(max(u));
            obj.movies=[];
        end
        
        function handler=scatterPlotSingleTimeStep(obj,K)
            if(obj.dim==2)
                handler=obj.scatterPlotSingleTimeStep2D(K);
            elseif(obj.dim==3)
                handler=obj.scatterPlotSingleTimeStep3D(K);
            else
                error('Only dimension 2 and 3 can be plotted.');
            end
        end
        
        function handler=scatterPlotSingleTimeStep2D(obj,K)
            if(not(obj.dim==2))
                error('The dimension of the mesh should be 2');
            end
            if (K<1 || K>length(obj.u(1,:)))
                error(['K has to be between 1 and ',length(obj.u(1,:))]);
            end
            handler=scatter(obj.mesh.points(:,1),obj.mesh.points(:,2),10.*ones(length(obj.mesh.points(:,1)),1),obj.colorConverter(K),'filled');
            obj.u(:,K)
        end
        
        function handler=scatterPlotSingleTimeStep3D(obj,K)
            if(not(obj.dim==3))
                error('The dimension of the mesh should be 3');
            end
            if (K<1 || K>length(obj.u(1,:)))
                error(['K has to be between 1 and ',length(obj.u(1,:))]);
            end
            handler=scatter3(obj.mesh.points(:,1),obj.mesh.points(:,2),obj.mesh.points(:,3),10.*ones(length(obj.mesh.points(:,1)),1),obj.colorConverter(K),'filled');
        end
        
        function movie=animateScatterPlot(obj)
            figure('Name','Scatter plot animation');
            obj.scatterPlotSingleTimeStep(1);
            axis tight manual;
            ax = gca;
            ax.NextPlot = 'replaceChildren';
            drawnow;
            movie(1) = getframe;
            [~,loops]=size(obj.u);
            movie(loops) = struct('cdata',[],'colormap',[]);
            for j = 2:loops
                obj.scatterPlotSingleTimeStep(j);
                drawnow;
                movie(j) = getframe;
            end
            obj.appendMovie(movie);
        end
        
        function movie=shrinkingPlot(obj,threshold)
        %Plots an animation over the timesteps in u where u<=threshold. The
        %animation will be the bounary of the points where u<=thresholds
        %holds.
            figure('Name','Shrinking plot animation');
            obj.shrinkingPlotSingleStep(threshold,1);
            axis tight manual;
            ax = gca;
            ax.NextPlot = 'replaceChildren';
            drawnow;
            movie(1) = getframe;
            [~,loops]=size(obj.u);
            movie(loops) = struct('cdata',[],'colormap',[]);
            for j = 2:loops
                obj.shrinkingPlotSingleStep(threshold,j);
                drawnow;
                movie(j) = getframe;
            end
            obj.appendMovie(movie);
        end
            
        
        function handler=shrinkingPlotSingleStep(obj,threshold,K)
        %Plots the boundary of the area/volume in which u<=threshold in
        %timestep k
            %Find the position of the values of u where u<=threshold
            position=find(obj.u(:,K)<=threshold);
            %If all points fullfill u<=threshold plot empty plot
            if(isempty(position))
                handler=obj.boundaryPlot([]);
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
                return
            end
            tri=triangulation(triangles,points);
            %Find boundry of the triangulation and plot it
            handler=obj.boundaryPlot(tri.freeBoundary());
        end
        
        function handler=boundaryPlot(obj,boundary)
            if(obj.dim==2)
                handler=obj.boundaryPlot2D(boundary);
            elseif(obj.dim==3)
                handler=obj.boundaryPlot3D(boundary);
            else
                error('Only dimension 2 and 3 can be plotted.');
            end
        end
        
        function handler=boundaryPlot2D(obj,boundary)
        %plots the given boundary within the 2D-mesh in obj. The boundary does
        %not need to be the boundary of the mesh.
            %Find path through the points that represent boundary
            if(isempty(boundary))
                %Plot nothing
                handler=scatter3(obj.mesh.points(1,1),obj.mesh.points(1,2),obj.mesh.points(1,3),'MarkerEdgeColor','white');
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
        %plots the given boundary within the 3D-mesh in obj. The boundary does
        %not need to be the boundary of the mesh.
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
        
        function [plain_points,projection]=getPlainPoints(obj,eps,E)
            projection=zeros(length(obj.mesh.points(:,1)),1);
            plain_points=zeros(length(obj.mesh.points(:,1)),3);
            counter=0;
            for k=1:length(obj.mesh.points(:,1))
                if(E.closeTo(obj.mesh.points(k,:)',eps))
                    counter=counter+1;
                    projection(counter)=k;
                    plain_points(counter,:)=obj.mesh.points(k,:);
                end
            end
            plain_points=plain_points(1:counter,:);
            projection=projection(1:counter);
        end 
        
        function movie=animateSlicePlot(obj,eps,varargin)
        %Plots an animation over the timesteps in u on the slice where x=0.
            if(nargin==2)
                E=Plain(eye(3),zeros(3,1));
            elseif(nargin==4)
                E=Plain(varargin{1},varargin{2});
            else
                error('animateSlicePlot needs 0 or 2 input values (next to the plotter object)');
            end
            [plain_points,projection]=getRotatedPlainPoints(obj,eps,E);
            
            figure('Name','Slice plot animation');
            obj.slicePlotSingleStep(plain_points,projection,1);
            axis tight manual;
            ax = gca;
            ax.NextPlot = 'replaceChildren';
            drawnow;
            movie(1) = getframe;
            [~,loops]=size(obj.u);
            movie(loops) = struct('cdata',[],'colormap',[]);
            for j = 2:loops
                obj.slicePlotSingleStep(plain_points,projection,j);
                drawnow;
                movie(j) = getframe;
            end
            obj.appendMovie(movie);
        end
        
        function [plain_points,projection]=getRotatedPlainPoints2(obj,eps,E)
            projection=zeros(length(obj.mesh.points(:,1)),1);
            plain_points=zeros(length(obj.mesh.points(:,1)),3);
            counter=0;
            for k=1:length(obj.mesh.points(:,1))
                if(E.closeTo(obj.mesh.points(k,:)',eps))
                    counter=counter+1;
                    projection(counter)=k;
                    plain_points(counter,:)=E.plain_inverse(obj.mesh.points(k,:)')';
                end
            end
            plain_points=plain_points(1:counter,2:3);
            projection=projection(1:counter);
        end 
        
        function handler=slicePlotSingleStep2(obj,plain_points,projection,K)
            C=Plotter.colorConverterStatic(obj.u(projection,K),obj.u_min,obj.u_max);
            handler=scatter(plain_points(:,1),plain_points(:,2),10.*ones(length(plain_points(:,1)),1),C,'filled');
        end
        
        function movie=animateSlicePlot2(obj,eps,varargin)
        %Plots an animation over the timesteps in u on the slice where x=0.
            if(nargin==2)
                E=Plain(eye(3),zeros(3,1));
            elseif(nargin==4)
                E=Plain(varargin{1},varargin{2});
            else
                error('animateSlicePlot needs 0 or 2 input values (next to the plotter object)');
            end
            [plain_points,projection]=getRotatedPlainPoints2(obj,eps,E);
            
            figure('Name','Slice plot animation');
            obj.slicePlotSingleStep2(plain_points,projection,1);
            axis tight manual;
            ax = gca;
            ax.NextPlot = 'replaceChildren';
            drawnow;
            movie(1) = getframe;
            [~,loops]=size(obj.u);
            movie(loops) = struct('cdata',[],'colormap',[]);
            for j = 2:loops
                obj.slicePlotSingleStep2(plain_points,projection,j);
                drawnow;
                movie(j) = getframe;
            end
            obj.appendMovie(movie);
        end
    end
end