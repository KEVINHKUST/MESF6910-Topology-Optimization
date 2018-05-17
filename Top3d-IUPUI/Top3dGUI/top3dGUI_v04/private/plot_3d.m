% DISPLAY 3D TOPOLOGY (ISO-VIEW)
% Special thanks to Tomas Zegard (UIUC) for quick plot function
function [p] = plot_3d(x, cutoff, Alpha, type, fcolor, disp)

if nargin < 6
    disp = x;
    if nargin < 5
        fcolor = 'r';
        if nargin < 4
            type = 'patch';
        end
    end
end

[nely,nelx,nelz] = size(x);

switch type
    case 'patch'
        x_logic = x>cutoff;
        
        aux = zeros(nely,nelx,nelz+1);
        aux(:,:,1) = x_logic(:,:,1); aux(:,:,nelz+1) = -x_logic(:,:,nelz);
        aux(:,:,2:nelz) = diff(x_logic,[],3);
        ind = find(aux);
        faceXY = zeros(numel(ind),4); faceXYcolor = zeros(numel(ind),1);
        [I,J,K] = ind2sub([nely,nelx,nelz+1],ind);
        faceXY(:,1) = sub2ind([nely+1 nelx+1 nelz+1],I,J,K);
        faceXY(:,2:end) = [faceXY(:,1)+nely+1 faceXY(:,1)+nely+2 faceXY(:,1)+1];
        ind_up = find(aux(ind)>0); ind_lo = find(aux(ind)<0);
        faceXYcolor(ind_up) = x(sub2ind([nely,nelx,nelz],I(ind_up),J(ind_up),K(ind_up)));
        faceXYcolor(ind_lo) = x(sub2ind([nely,nelx,nelz],I(ind_lo),J(ind_lo),K(ind_lo)-1));
        
        aux = zeros(nely,nelx+1,nelz);
        aux(:,1,:) = x_logic(:,1,:); aux(:,nelx+1,:) = -x_logic(:,nelx,:);
        aux(:,2:nelx,:) = diff(x_logic,[],2);
        ind = find(aux);
        faceYZ = zeros(numel(ind),4); faceYZcolor = zeros(numel(ind),1);
        [I,J,K] = ind2sub([nely,nelx+1,nelz],ind);
        faceYZ(:,1) = sub2ind([nely+1 nelx+1 nelz+1],I,J,K);
        faceYZ(:,2:end) = [faceYZ(:,1)+1 faceYZ(:,1)+(nelx+1)*(nely+1)+1 faceYZ(:,1)+(nelx+1)*(nely+1)];
        ind_up = find(aux(ind)>0); ind_lo = find(aux(ind)<0);
        faceYZcolor(ind_up) = x(sub2ind([nely,nelx,nelz],I(ind_up),J(ind_up),K(ind_up)));
        faceYZcolor(ind_lo) = x(sub2ind([nely,nelx,nelz],I(ind_lo),J(ind_lo)-1,K(ind_lo)));
        
        aux = zeros(nely+1,nelx,nelz);
        aux(1,:,:) = x_logic(1,:,:); aux(nely+1,:,:) = -x_logic(nely,:,:);
        aux(2:nely,:,:) = diff(x_logic,[],1);
        ind = find(aux);
        faceZX = zeros(numel(ind),4); faceZXcolor = zeros(numel(ind),1);
        [I,J,K] = ind2sub([nely+1,nelx,nelz],ind);
        faceZX(:,1) = sub2ind([nely+1 nelx+1 nelz+1],I,J,K);
        faceZX(:,2:end) = [faceZX(:,1)+(nelx+1)*(nely+1) faceZX(:,1)+(nelx+2)*(nely+1) faceZX(:,1)+(nely+1)];
        ind_up = find(aux(ind)>0); ind_lo = find(aux(ind)<0);
        faceZXcolor(ind_up) = x(sub2ind([nely,nelx,nelz],I(ind_up),J(ind_up),K(ind_up)));
        faceZXcolor(ind_lo) = x(sub2ind([nely,nelx,nelz],I(ind_lo)-1,J(ind_lo),K(ind_lo)));
        
        cla, hold on, view(30,30), rotate3d on, axis equal, axis([0 nelx 0 nelz 0 nely]), box
        set(gca,'YDir','reverse','ZDir','reverse','ZtickLabel',flipud(get(gca,'Ztick')'));
        % set(gcf,'MenuBar','none','ToolBar','none')
        [X,Y,Z] = meshgrid(0:nelx,0:nely,0:nelz);
        p = patch('Faces',[faceXY; faceYZ; faceZX],'Vertices',[X(:) Z(:) Y(:)],...
            'FaceColor','flat','FaceVertexCData',[faceXYcolor; faceYZcolor; faceZXcolor],'FaceAlpha',Alpha);
        colormap(flipud(gray)), caxis([0 1]), colorbar('off'), drawnow
    case 'isosurface'
        aux = zeros(nely+2,nelx+2,nelz+2);
        aux(2:end-1,2:end-1,2:end-1) = x;

        cla, hold on, view(30,30), rotate3d on, axis equal, axis([0 nelx 0 nelz 0 nely]), box
        set(gca,'YDir','reverse','ZDir','reverse','ZtickLabel',flipud(get(gca,'Ztick')'));
        % set(gcf,'MenuBar','none','ToolBar','none')
        [X,Y,Z] = meshgrid(0:nelx+1,0:nely+1,0:nelz+1);
        p = patch(isosurface(X-0.5,Z-0.5,Y-0.5,aux,cutoff),...
            'FaceColor',fcolor,'EdgeColor','none','FaceAlpha',Alpha);
        camlight, lighting gouraud;
        colorbar('off'), drawnow
    case 'isonormals'
        aux = zeros(nely+2,nelx+2,nelz+2);
        aux(2:end-1,2:end-1,2:end-1) = x;

        cla, hold on, view(30,30), rotate3d on, axis equal, axis([0 nelx 0 nelz 0 nely]), box
        set(gca,'YDir','reverse','ZDir','reverse','ZtickLabel',flipud(get(gca,'Ztick')'));
        % set(gcf,'MenuBar','none','ToolBar','none')
        [X,Y,Z] = meshgrid(0:nelx+1,0:nely+1,0:nelz+1);
        p = patch(isosurface(X-0.5,Z-0.5,Y-0.5,aux,cutoff),...
            'FaceColor',fcolor,'EdgeColor','none','FaceAlpha',Alpha);
        isonormals(X,Y,Z,aux,p)
        camlight, lighting gouraud;
        colorbar('off'), drawnow
    case 'disp'
        
        x_logic = x>cutoff;
        
        aux = zeros(nely,nelx,nelz+1);
        aux(:,:,1) = x_logic(:,:,1); aux(:,:,nelz+1) = -x_logic(:,:,nelz);
        aux(:,:,2:nelz) = diff(x_logic,[],3);
        ind = find(aux);
        faceXY = zeros(numel(ind),4); faceXYcolor = zeros(numel(ind),1);
        [I,J,K] = ind2sub([nely,nelx,nelz+1],ind);
        faceXY(:,1) = sub2ind([nely+1 nelx+1 nelz+1],I,J,K);
        faceXY(:,2:end) = [faceXY(:,1)+nely+1 faceXY(:,1)+nely+2 faceXY(:,1)+1];
        ind_up = find(aux(ind)>0); ind_lo = find(aux(ind)<0);
        faceXYcolor(ind_up) = disp(sub2ind([nely,nelx,nelz],I(ind_up),J(ind_up),K(ind_up)));
        faceXYcolor(ind_lo) = disp(sub2ind([nely,nelx,nelz],I(ind_lo),J(ind_lo),K(ind_lo)-1));
        
        aux = zeros(nely,nelx+1,nelz);
        aux(:,1,:) = x_logic(:,1,:); aux(:,nelx+1,:) = -x_logic(:,nelx,:);
        aux(:,2:nelx,:) = diff(x_logic,[],2);
        ind = find(aux);
        faceYZ = zeros(numel(ind),4); faceYZcolor = zeros(numel(ind),1);
        [I,J,K] = ind2sub([nely,nelx+1,nelz],ind);
        faceYZ(:,1) = sub2ind([nely+1 nelx+1 nelz+1],I,J,K);
        faceYZ(:,2:end) = [faceYZ(:,1)+1 faceYZ(:,1)+(nelx+1)*(nely+1)+1 faceYZ(:,1)+(nelx+1)*(nely+1)];
        ind_up = find(aux(ind)>0); ind_lo = find(aux(ind)<0);
        faceYZcolor(ind_up) = disp(sub2ind([nely,nelx,nelz],I(ind_up),J(ind_up),K(ind_up)));
        faceYZcolor(ind_lo) = disp(sub2ind([nely,nelx,nelz],I(ind_lo),J(ind_lo)-1,K(ind_lo)));
        
        aux = zeros(nely+1,nelx,nelz);
        aux(1,:,:) = x_logic(1,:,:); aux(nely+1,:,:) = -x_logic(nely,:,:);
        aux(2:nely,:,:) = diff(x_logic,[],1);
        ind = find(aux);
        faceZX = zeros(numel(ind),4); faceZXcolor = zeros(numel(ind),1);
        [I,J,K] = ind2sub([nely+1,nelx,nelz],ind);
        faceZX(:,1) = sub2ind([nely+1 nelx+1 nelz+1],I,J,K);
        faceZX(:,2:end) = [faceZX(:,1)+(nelx+1)*(nely+1) faceZX(:,1)+(nelx+2)*(nely+1) faceZX(:,1)+(nely+1)];
        ind_up = find(aux(ind)>0); ind_lo = find(aux(ind)<0);
        faceZXcolor(ind_up) = disp(sub2ind([nely,nelx,nelz],I(ind_up),J(ind_up),K(ind_up)));
        faceZXcolor(ind_lo) = disp(sub2ind([nely,nelx,nelz],I(ind_lo)-1,J(ind_lo),K(ind_lo)));
        
        cla, hold on, view(30,30), rotate3d on, axis equal, axis([0 nelx 0 nelz 0 nely]), box, colorbar('off')
        set(gca,'YDir','reverse','ZDir','reverse','ZtickLabel',flipud(get(gca,'Ztick')'));
        % set(gcf,'MenuBar','none','ToolBar','none')
        [X,Y,Z] = meshgrid(0:nelx,0:nely,0:nelz);
        p = patch('Faces',[faceXY; faceYZ; faceZX],'Vertices',[X(:) Z(:) Y(:)],...
            'FaceColor','flat','FaceVertexCData',[faceXYcolor; faceYZcolor; faceZXcolor],'FaceAlpha',Alpha);
        colormap(jet), caxis([min(disp(:)) max(disp(:))]), drawnow, %colorbar('location','southoutside')
end