function DwPlot(A,xyRes,xzRes,varargin)
%DWPLOT plot the inscribed approaximation of the DW shell. This function
%relies on the function DWSHELL, FOVALS, make sure you have them in the MATLAB path.
%You are required to specify two rotational resolution, i.e., number of sample
%points in each rotating direction. You can also specify the plotting
%specifications, such as 'EdgeColor', 'FaceAlpha'(degree of transparency),
%etc. as long as it is supported by the built-in trisurf function.
% DwPlotlot(A,100,100) plot the DW shell of A with angular resolution
% 100(xy) x 100 (xz)
% DwPlot(A,100,100,'FaceAlpha',0.5) adjust the face transparency to 0.5
% DwPlot(A,100,120,'EdgeColor','none') hide the face edges

%% input processing
    plotSpec = {'FaceAlpha',0.7}; % default face transparency
    if nargin > 3
         plotSpec = varargin;
    end
    hold on;

%% plotting the shell
    bp = DwShell(A,xyRes,xzRes);
    x = bp(1,:); y = bp(2,:); z = bp(3,:);
    [k1,av1] = convhull(x,y,z);
    trisurf(k1,x,y,z,plotSpec{:});

%% find the bounds on the DW shell
    % find the rectangular which bounds the spectrum of the matrix
    HA = (A+A')/2; SA = (A-A')/(2*1i);
    eigHA = eig(HA); eigSA = eig(SA);
    xMax = max(eigHA); xMin = min(eigHA);
    yMax = max(eigSA); yMin = min(eigSA); % Bendixson's rectangular bound on numerical range
    
    % plot the parabola over a larger rectangle area
    xRelax = (xMax-xMin)/4;
    yRelax = (yMax-yMin)/4;  
    % parabola bound
    xRes = 500; yRes = 500; nBound = 1000;  
    [X,Y] = meshgrid(linspace(xMin-xRelax,xMax+xRelax,xRes),...
        linspace(yMin-yRelax,yMax+yRelax,yRes));
    Z = X.^2 + Y.^2;
    surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.5);

%     % lid of the pillar
    sigma1 = max(svd(A));
    bpFoval = fovals(A,nBound); % boundary point of numerical range
    lid = poly2mask((bpFoval(1,:)-xMin+xRelax)/(X(1,2)-X(1,1)),...
        (bpFoval(2,:)-yMin+yRelax)/(Y(2,1)-Y(1,1)),xRes,yRes);
    lidZ = sigma1^2*lid;
    surf(X,Y,lidZ,'EdgeColor','none','FaceAlpha',.2);
    
%% post process
    xlabel('$\Re (x^*Ax)$','Interpreter','latex');
    ylabel('$\Im (x^*Ax)$','Interpreter','latex');
    zlabel('$\| Ax \|^2$','Interpreter','latex');
    allAxes = findall(0,'type','axes');
    for i = 1:length(allAxes)
        box(allAxes(i), 'on');
        allAxes(i).FontSize = 20;
    end
end

