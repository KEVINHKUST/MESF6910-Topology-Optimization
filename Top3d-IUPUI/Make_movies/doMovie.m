function doMovie
% Make a topology optimizatin movie
% By K. Liu
 
% cutoffL allow smooth animation from 1 to 50 iterations.
% change 0.5 to the desired element density cutoff value
% change smaller 50 to allow faster animation, larger 50 to allow slower animation
cutoffL = linspace(0,0.5,50);
V = VideoWriter('myTopopt.avi','Uncompressed AVI');
open(V);
% Do Movie 
i = 0;
while (1)
    i = i + 1;
    filename = ['xPhys_it',num2str(i,'%03d'),'.mat'];
    if exist(filename,'file') == 0
        break;
    end
    load(filename)
    if i <= length(cutoffL)
        cutoff = cutoffL(i);
    else
        cutoff = max(cutoffL);
    end
     
    plot_3d(xPhys, cutoff, 1, [0.7 0.7 0.7]);
    % Syntax
    % plot_3d
    % [p] = plot_3d(x, cutoff, Alpha, fcolor)
    % p = patch(isosurface(X-0.5,Z-0.5,Y-0.5,aux,cutoff) 'FaceColor',fcolor,
    % 'EdgeColor','none','FaceAlpha',Alpha);
    
    drawnow
    mov(i) = getframe(gcf); %#ok<AGROW>
    % It suppresses mlint warnings. In this specific case, it is about not 
    % pre-allocating an array.
    % mlint is one of the static code analysis tools that Matlab has. It 
    % finds possible errors and shows warnings.
    
    % mlint is replaced with checkcode
    % checkcode:Check MATLAB code files for possible problems
    % checkcode(filename) displays messages about filename that report potential 
    % problems and opportunities for code improvement. These messages are 
    % sometimes referred to as Code Analyzer messages. The line number in 
    % the message is a hyperlink that you can click to go directly to that 
    % line in the Editor. The exact text of the checkcode messages is subject 
    % to some change between versions.
end
writeVideo(V,mov);
% WriteVideo(v,frame) writes one or more movie frames typically returned by
% the getframe function.
close(V);

end
 
% DISPLAY 3D TOPOLOGY (ISO-VIEW)
% Special thanks to Tomas Zegard (UIUC) for quick plot function
% Code snippets from Top3d (top3dapp.com)
% Top3d - An efficient 3D topology optimization program written in MATLAB
% By K. Liu and A. Tovar
% http://www.top3dapp.com
function [p] = plot_3d(x, cutoff, Alpha, fcolor)
 
if nargin < 4, fcolor = 'r'; end
% nargin returns the number of function input arguments given in the call 
% to the currently executing function. Use this syntax in the body of a 
% function only.
 
[nely,nelx,nelz] = size(x);
aux = zeros(nely+2,nelx+2,nelz+2);
aux(2:end-1,2:end-1,2:end-1) = x;
 
cla, hold on, view(30,30), rotate3d on, axis equal, axis([0 nelx 0 nelz 0 nely]), box
% 1.clear axes
% 2.view(az,el) and view([az,el]) set the viewing angle for a three-dimensional 
% plot. 
% The azimuth, az, is the horizontal rotation about the z-axis as measured 
% in degrees from the negative y-axis. Positive values indicate 
% counterclockwise rotation of the viewpoint. 
% The elevation, el is the vertical elevation of the viewpoint in degrees. Positive values 
% of elevation correspond to moving above the object; negative values correspond to moving below the object.
% 3.rotate3d on:enables mouse-base rotation on all axes within the current figure.
% 4.axis equal: Use axis lines with equal lengths. Adjust the increments between data units accordingly.
% 5.Set the limits for the x,y,z-axis 
% 6.box:toggles the display of the box outline.

set(gca,'YDir','reverse','ZDir','reverse','ZtickLabel',flipud(get(gca,'Ztick')'));
% set(gcf,'MenuBar','none','ToolBar','none')
% Set graphics object properties
[X,Y,Z] = meshgrid(0:nelx+1,0:nely+1,0:nelz+1);
p = patch(isosurface(X-0.5,Z-0.5,Y-0.5,aux,cutoff),...
    'FaceColor',fcolor,'EdgeColor','none','FaceAlpha',Alpha);
camlight, lighting gouraud;
% camlight with no arguments is the same as camlight('right').
% camlight('right') creates a light right and up from camera.
% lighting gouraud calculates the vertex normals and interpolates 
% linearly across the faces. Select this method to view curved surfaces.




axis off; box off; set(gcf, 'color', [1 1 1])
drawnow
end