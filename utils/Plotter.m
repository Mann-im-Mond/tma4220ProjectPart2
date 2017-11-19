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
            cC=@(v) [to01Convert(v),0.*v,1-to01Convert(v)];
            rgb=cC(u);
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
            if not(N==length(mesh.points))
                error('The number of rows must be equal to the number of points in mesh');
            end
            if (M<1)
                error('u must not be empty');
            elseif (M==1)
                timeInterval.t_max=timeInterval.t_max_0;
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
            if (K<1 || K>length(obj.u'))
                error(['K has to be between 1 and ',length(obj.u')]);
            end
            handler=scatter(obj.mesh.points(:,1),obj.mesh.points(:,2),10.*ones(length(obj.mesh.points),1),obj.colorConverter(K),'filled');
            obj.u(:,K)
        end
        
        function handler=scatterPlotSingleTimeStep3D(obj,K)
            if(not(obj.dim==3))
                error('The dimension of the mesh should be 3');
            end
            if (K<1 || K>length(obj.u'))
                error(['K has to be between 1 and ',length(obj.u')]);
            end
            handler=scatter3(obj.mesh.points(:,1),obj.mesh.points(:,2),obj.mesh.points(:,3),10.*ones(length(obj.mesh.points),1),obj.colorConverter(K),'filled');
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
                obj.shrinkingPlotSingleStep(threshold,k);
                drawnow;
                movie(j) = getframe;
            end
            obj.appendMovie(movie);
        end
            
        
        function handler=shrinkingPlotSingleStep(obj,threshold,k)
            position=find(obj.u(:,k)<=threshold);
            projection=zeros(length(obj.u(:,k)),1);
            counter=1;
            for k=1:length(obj.u(:,k))
                if(position(counter)==k)
                    projection(k)=counter;
                    counter=counter+1;
                end
            end
            points=obj.mesh.points(position,:);
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
            tri=triangulation(triangles,points);
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
            faces=zeros(1,length(boundary));
            faces(1)=boundary(1,1);
            positionX=1;
            positionY=1;
            for k=2:length(boundary)
                faces(k)=boundary(positionY,3-positionX);
                next=find(boundary(:,positionX)==faces(k));
                if(isempty(next))
                    positionX=3-positionX;
                    next=find(boundary(:,positionX)==faces(k));
                    positionY=next(next~=positionY);
                elseif(length(next)==1)
                    positionY=next;
                else
                    error('The boundary is not a circle.');
                end
            end
            if(faces(1)~=boundary(positionY,3-positionX))
                error('The boundary is not a circle.');
            end
            handler=patch('Faces',faces,'Vertices',obj.mesh.points,'EdgeColor','blue','FaceColor','none');
        end
        
        function handler=boundaryPlot3D(obj,boundary)
            l=length(boundary);
            X=zeros(3,l);
            Y=zeros(3,l);
            Z=zeros(3,l);
            C=zeros(3,l);
            C(:,3)=1;
            for k=1:l
                X(:,k)=obj.mesh.points(boundary(k,:),1);
                Y(:,k)=obj.mesh.points(boundary(k,:),2);
                Z(:,k)=obj.mesh.points(boundary(k,:),3);
            end
            handler=fill3(X,Y,Z,C);
        end
    end
end