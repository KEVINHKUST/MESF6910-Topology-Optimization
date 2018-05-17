% A 169 LINE 3D TOPOLOGY OPITMIZATION CODE BY LIU AND TOVAR (JUL 2013)
function [p, xPhys, disp] = top3d(handles)
% DEMO
volfrac = handles.domain.volfrac;
rmin = handles.optpar.rmin;
nelx = handles.domain.nelx;
nely = handles.domain.nely;
nelz = handles.domain.nelz;
% USER-DEFINED LOOP PARAMETERS
maxloop = handles.optpar.maxloop;    % Maximum number of iterations
tolx = handles.optpar.tolx;      % Termination criterion
displayflag = handles.optpar.displayflag;  % Display structure flag
cutoff  = 0.5;
% USER-DEFINED MATERIAL PROPERTIES
E0    = handles.matlcon.E0;    % Young's modulus of solid material
Emin  = 1e-9;                  % Young's modulus of void-like material
nu    = handles.matlcon.nu;    % Poisson's ratio
penal = handles.matlcon.p;     % Penalization Power
% USER-DEFINED LOAD DOFs
% il = nelx; jl = 0; kl = 0:nelz;                         % Coordinates
% loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
% loaddof = 3*loadnid(:) - 1;                             % DOFs
loaddof = handles.lc.loaddof;
loadscale = handles.lc.loadscale;
% USER-DEFINED SUPPORT FIXED DOFs
% [jf,kf] = meshgrid(1:nely+1,1:nelz+1);                  % Coordinates
% fixednid = (kf-1)*(nely+1)*(nelx+1)+jf;                 % Node IDs
fixeddof = handles.bc.fixeddof; % DOFs
% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,loadscale,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
KE = lk_H8(nu);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = kron(edofMat,ones(24,1))';
jK = kron(edofMat,ones(1,24))';
% PREPARE FILTER
step = ceil(rmin)-1;
iH = zeros(nele*(2*step+1)^3,1);
jH = zeros(size(iH)); vH = zeros(size(iH));
n = 0;
for el=1:nele
    [i,j,k] = ind2sub([nely,nelx,nelz],el);
    [ispan,jspan,kspan] = meshgrid(max(1,i-step):min(nely,i+step),max(1,j-step):min(nelx,j+step),max(1,k-step):min(nelz,k+step));
    dist = max(0,rmin-sqrt((ispan-i).^2 + (jspan-j).^2 + (kspan-k).^2));
    vH(n+(1:numel(dist))) = dist(:);
    iH(n+(1:numel(dist))) = el;
    jH(n+(1:numel(dist))) = sub2ind([nely nelx nelz],ispan,jspan,kspan);
    n = n + numel(dist);
end
iH(n+1:end)=[]; jH(n+1:end)=[]; vH(n+1:end)=[];
H = sparse(iH,jH,vH);
Hs = sum(H,2);
% INITIALIZE ITERATION
x = repmat(volfrac,[nely,nelx,nelz]);
xPhys = x;
loop = 0;
change = 1;
% START ITERATION
while change > tolx && loop < maxloop
    loop = loop+1;
    % FE-ANALYSIS
    sK = KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin));
    K = sparse(iK(:),jK(:),sK(:)); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    c = sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)));
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx,nelz);
    % FILTERING AND MODIFICATION OF SENSITIVITIES
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
    % OPTIMALITY CRITERIA UPDATE
    if loop<15, l1 = 0.0;      l2 = 1e9;      move = 0.15;
    else        l1 = lmid/1.1; l2 = lmid*1.1; move = 0.10;
    end
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
        xPhys(:) = (H*xnew(:))./Hs;
        if sum(xPhys(:)) > volfrac*nele, l1 = lmid; else l2 = lmid; end
    end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    
    str = sprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
    title(str); fprintf(1, str);
    % PRINT RESULTS
    if displayflag
        % PLOT DENSITIES
        cla, plot_3d(xPhys,0.5,1);
    end
end
set(handles.axes_res,'Visible','on');
cla; p = plot_3d(xPhys,cutoff,1);
% Maniuplate displacements
disp = zeros(nele, 3);
for i = 1:nele
    disp(i,:) =  sum(reshape(U(edofMat(i, :)),3,[]),2)';
end
% Resulant displacement
disp = reshape(sum(disp.^2, 2), nely, nelx, nelz); % resulant displacement
end
% ===================== AUXILIARY FUNCTIONS ===============================
% GENERATE ELEMENT STIFFNESS MATRIX
function [KE] = lk_H8(nu)
A = [ 32   6  -8   6  -6   4   3  -6 -10   3  -3  -3  -4  -8;
     -48   0   0 -24  24   0   0   0  12 -12   0  12  12  12];
k = (1/144)*A'*[1; nu];
% GENERATE SIX SUB-MATRICES AND THEN GET KE MATRIX
K1 = k([ 1  2  2  3  5  5;
         2  1  2  4  6  7;
         2  2  1  4  7  6;
         3  4  4  1  8  8;
         5  6  7  8  1  2;
         5  7  6  8  2  1]);
K2 = k([ 9  8 12  6  4  7;
         8  9 12  5  3  5;
        10 10 13  7  4  6;
         6  5 11  9  2 10;
         4  3  5  2  9 12;
        11  4  6 12 10 13]);
K3 = k([ 6  7  4  9 12  8;
         7  6  4 10 13 10;
         5  5  3  8 12  9;
         9 10  2  6 11  5;
        12 13 10 11  6  4;
         2 12  9  4  5  3]);
K4 = k([14 11 11 13 10 10;
        11 14 11 12  9  8;
        11 11 14 12  8  9;
        13 12 12 14  7  7;
        10  9  8  7 14 11;
        10  8  9  7 11 14]);
K5 = k([ 1  2  8  3  5  4;
         2  1  8  4  6 11;
         8  8  1  5 11  6;
         3  4  5  1  8  2;
         5  6 11  8  1  8;
         4 11  6  2  8  1]);
K6 = k([14 11  7 13 10 12;
        11 14  7 12  9  2;
         7  7 14 10  2  9;
        13 12 10 14  7 11;
        10  9  2  7 14  7;
        12  2  9 11  7 14]);
KE = 1/((nu+1)*(1-2*nu))*[ K1  K2  K3  K4 ;
                           K2' K5  K6  K3';
                           K3' K6  K5' K2';
                           K4  K3  K2  K1'];
end

% =========================================================================
% === This code was written by K Liu and A Tovar, Dept. of Mechanical   ===
% === Engineering, Indiana University-Purdue University Indianapolis,   ===
% === Indiana, United States of America                                 ===
% === ----------------------------------------------------------------- ===
% === Please send your suggestions and comments to: kailiu@iupui.edu    ===
% === ----------------------------------------------------------------- ===
% === The code is intended for educational purposes, and the details    ===
% === and extensions can be found in the paper:                         ===
% === '' ''                                                             ===
% === ----------------------------------------------------------------- ===
% === The code as well as an uncorrected version of the paper can be    ===
% === downloaded from the website: http://engr.iupui.edu/~tovara/top3d  ===
% === ----------------------------------------------------------------- ===
% === Disclaimer:                                                       ===
% === The authors reserves all rights for the program.                  ===
% === The code may be distributed and used for educational purposes.    ===
% === The authors do not guarantee that the code is free from errors, and =
% === they shall not be liable in any event caused by the use of the code.=
% =========================================================================