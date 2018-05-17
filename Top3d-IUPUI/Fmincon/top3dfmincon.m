% AN 169 LINE 3D TOPOLOGY OPITMIZATION CODE BY LIU AND TOVAR (JUL 2013)
function top3dfmincon(nelx,nely,nelz,volfrac,penal,rmin)
% nelx,nely,nelz,volfrac,penal span multiple functions See
% reference:https://ww2.mathworks.cn/help/matlab/matlab_prog/nested-functions.html?lang=en

% Variables within nested functions are accessible to more than just their
% immediate function. A variable, x, to which you assign a value or use
% within a nested function resides in the workspace of the outermost
% function that both contains the nested function and accesses x.
% Therefore, the scope of x is the function to which this workspace
% belongs, and all functions nested to any level within that function.

% top3dfmincon is the parent function
% myObjFcn,myConstrFcn,myHessianFcn,myOutputFcn are the nested functions

% USER-DEFINED LOOP PARAMETERS
maxloop = 200;    % Maximum number of iterations
tolx = 0.01;      % Terminarion criterion
displayflag = 0;  % Display structure flag
% USER-DEFINED MATERIAL PROPERTIES
E0 = 1;           % Young's modulus of solid material
Emin = 1e-9;      % Young's modulus of void-like materiald
nu = 0.3;         % Poisson's ratio
% USER-DEFINED LOAD DOFs
[il,jl,kl] = meshgrid(nelx, 0, 0:nelz);                 % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
% USER-DEFINED SUPPORT FIXED DOFs
[iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-1,ndof,1);
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
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
% PREPARE FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
% INITIALIZE ITERATION
x = repmat(volfrac,[nely,nelx,nelz]);
xPhys = x; 
%% Step.1: Initialize Fmincon
global ce 

% ce will be used in both myObjFcn and myHessianFcn functions 
% global var1... varN declares variables var1 ... varN as global in scope. 
% Ordinarily, each MATLAB? function has its own local variables, which are separate
% from those of other functions and from those of the base workspace.
% However, if several functions all declare a particular variable name as
% global, then they all share a single copy of that variable. Any change of
% value to that variable, in any function, is visible to all the functions
% that declare it as global. 
% If a variable with the same name as the global variable already exists in
% the current workspace, MATLAB issues a warning and changes the value of 
% that variable and its scope to match the global variable.

A = [];
B = [];
Aeq = [];
Beq = [];
LB = zeros(size(x));
UB = ones(size(x));
OPTIONS = optimset('TolX',tolx,...
    'MaxIter',maxloop,...
    'Algorithm','interior-point',...
    'GradObj','on', ...
    'GradConstr','on',...
    'Hessian','user-supplied',...
    'HessFcn',@myHessianFcn,...
    'Display','none',...
    'OutputFcn',@(x,optimValues,state) myOutputFcn(x,optimValues,state,displayflag),...
    'PlotFcns',@optimplotfval);

% optimset: Create or edit optimization options structure
% Optimset Reference: Optimization Options Reference/ Current and Legacy 
% Option Name Tables

% The function optimset creates an options structure that you can pass as 
% an input argument to the following four MATLAB? optimization
% functions:1.fminbnd:Find minimum of single-variable function on fixed 
%             interval 
%           2.fminsearch:Find minimum of unconstrained multivariable function
%            using derivative-free method
%           3.fzero:Root of nonlinear function 
%           4.lsqnonneg:Solve nonnegative linear least-squares problem 
% Alogrithm: "interior-point" chosen from "sqp","active-set","trust-region-reflective"
% Hessian: user-supplied Hessian of the Lagrange function
% GradObj, GradConstr, Hessian, HessFcn: With analytic expressions of the
% gradients of the objective function, constraints and Hessian of the
% Lagrange function, the optimization algorithm will execute faster than
% numerical approximation. 
% Display: display option is turned off. The iteration information will be 
% displaced in the OutputFcn. 
% OutputFcn: Output function is defined in myOutputFcn, which will display 
% iteration information and structural topology corresponding. 
% PlotFcns: Matlab provides some built-in plot function, @optimplotfval 
% plots the function value.

%% Step.2: Define Objective function
function [f, gradf] = myObjFcn(x) % The function to be minimized
    xPhys(:) = (H*x(:))./Hs;
    % FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),24*24*nele,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    c = sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)));
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    % FILTERING AND MODIFICATION OF SENSITIVITIES
    dc(:) = H*(dc(:)./Hs);
    % RETURN
    f = c; % Objective function:compliance
    gradf = dc(:); % Gradient function: dc 
end % myfun

%% Step.3: Define Hessian
function h = myHessianFcn(x, lambda)
% Reference: fmincon Interior-Point Algorithm with Analytic Hessian
% x:point x 
% lambda:Lagrange multiplier structure lambda
% 
    xPhys = reshape(x,nely,nelx,nelz);
    % Compute Hessian of Obj.
    Hessf = 2*(penal*(E0-Emin)*xPhys.^(penal-1)).^2 ./ (E0 + (E0-Emin)*xPhys.^penal) .* ce;
    Hessf(:) = H*(Hessf(:)./Hs);
    % Compute Hessian of constraints
    Hessc = 0; % Linear constraint
    % Hessian of Lagrange
    h = diag(Hessf(:)) + lambda.ineqnonlin*Hessc;
end % myHessianFcn

%% Step.4: Define Constraint
function [cneq, ceq, gradc, gradceq] = myConstrFcn(x)
    xPhys(:) = (H*x(:))./Hs;
    % Non-linear Constraints
    cneq  = sum(xPhys(:)) - volfrac*nele; % nonlinear inequality constraints c(x)<=0 
    gradc = ones(nele,1); % gradient of c(x)
    % Linear Constraints
    ceq     = []; % nonlinear equality constraints: ceq(x) = 0
    gradceq = []; % gradient of ceq(x)
end % mycon

%% Step.5: Define Output Function
function stop = myOutputFcn(x,optimValues,state,displayflag)
    stop = false;
    switch state
        case 'iter'
            % Make updates to plot or guis as needed
            xPhys = reshape(x, nely, nelx, nelz);
            %% PRINT RESULTS
            fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',optimValues.iteration,optimValues.fval, ...
                mean(xPhys(:)),optimValues.stepsize);
            %% PLOT DENSITIES
            if displayflag, figure(10); clf; display_3D(xPhys); end
            title([' It.:',sprintf('%5i',optimValues.iteration),...
                ' Obj. = ',sprintf('%11.4f',optimValues.fval),...
                ' ch.:',sprintf('%7.3f',optimValues.stepsize)]);
        case 'init'
            % Setup for plots or guis
            if displayflag
                figure(10)
            end
        case 'done'
            % Cleanup of plots, guis, or final plot
            figure(10); clf; display_3D(xPhys);
        otherwise
    end % switch
end % myOutputFcn

%% Step.6: Call Fmincon

fmincon(@myObjFcn, x, A, B, Aeq, Beq, LB, UB, @myConstrFcn, OPTIONS);

end

% === GENERATE ELEMENT STIFFNESS MATRIX ===
function [KE] = lk_H8(nu)
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/144*A'*[1; nu];

K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end
% === DISPLAY 3D TOPOLOGY (ISO-VIEW) ===
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.5)  % User-defined display density threshold
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
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
% === K. Liu and A. Tovar, "An efficient 3D topology optimization code  ===
% === written in Matlab", Struct Multidisc Optim, 50(6): 1175-1196, 2014, =
% === doi:10.1007/s00158-014-1107-x                                     ===
% === ----------------------------------------------------------------- ===
% === The code as well as an uncorrected version of the paper can be    ===
% === downloaded from the website: http://www.top3dapp.com/             ===
% === ----------------------------------------------------------------- ===
% === Disclaimer:                                                       ===
% === The authors reserves all rights for the program.                  ===
% === The code may be distributed and used for educational purposes.    ===
% === The authors do not guarantee that the code is free from errors, a