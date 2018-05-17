%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
%% Revised by ChenKaiwen HKUST
%  Added comments
% 26Mar2018

%% The main function
function top88(nelx,nely,volfrac,penal,rmin,ft)
%The filter radius rmin equals 0.04 times the width of the design domain in
%the paper
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]); % The element 
% Stiffness matrix k0 for an element with unit Young's modulus

nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
% node number index, each node has one number
% nodenrs matrix has (1+nely) rows, (1+nelx) columns
% The element is from 1 to (1+nelx)*(1+nely), which means there are
% (1+nelx)*(1+nely) nodes

% reshape function
% Syntax: B = reshape(A,sz1,...,szN)
% B = reshape(A,sz1,...,szN) reshapes A into a sz1-by-...-by-szN array where 
% sz1,...,szN indicates the size of each dimension.

edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
% Determine the first DOF index for all elements
% In total, there are nelx*nely elements
% the first DOF index means the x-direction DOF INDEX of the bottom-left node in
% every element

edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
% Assembly of the stiffness matrix
% The i-th row of edofMat contains the eight Degree of freedom indices
% corresponding to the i-th element, in terms of column-wise from top to
% bottom, from left to right
% Determine the eight DOF indices for each element

% repmat function: repeat copies of array
% Syntax:B = repmat(A,r1,...,rN)
% B = repmat(A,r1,...,rN) specifies a list of scalars,r1,...,rN, that
% describes how copies of  A are arranged in each dimension.
% Example:
% repmat([1 2; 3 4],2,3) = 
%      1     2     1     2     1     2
%      3     4     3     4     3     4
%      1     2     1     2     1     2
%      3     4     3     4     3     4

% repmat(edofVec,1,8) copies the coloumn vector into 8 coloumns which share
% the same number in each row,i.e. returns a matrix with eight columns
% which are all copies of the vector edofVec.
% Then the repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1) orders/gives 
% every element the right indices order based on the first DOF index.It 
% relates the indices of the eight DOFs of an element to the index of its 
% first DOF stored in the vector edofVec

iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
% iK: The row indices of the non-zero matrix entries
% iK is structured so that the indices iK(k) correspond to the (i,j)th
% entry of the stiffness matrix for element e, where k = i+8(j-1)+64(e-1)
% For every element, there are 8 degree of freedom. Each multiply the
% unit vector of length 8. 8*8=64. So in total there are 64*nelx*nely

% kron function: Kronecker tensor product
% Syntax: K = kron(A,B)
% K = kron(A,B) returns the Kronecker tensor product of matrices A and B. 
% If A is an m-by-n matrix and B is a p-by-q matrix, then kron(A,B) is an 
% m*p-by-n*q matrix formed by taking all possible products between the elements 
% of A and the matrix B.
% 

jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% jK: The column indices of the non-zero matrix entries
% jK(k) correspond to the (i, j )-th entry of the stiffness matrix for 
% element e, where k = i + 8( j ? 1) + 64(e ? 1).

% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union((1:2:2*(nely+1)),(2*(nelx+1)*(nely+1)));
alldofs = (1:2*(nely+1)*(nelx+1));
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
% Preallocated iH,jH and sH
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1; % 
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
% Hs and H remain invariant during the optimization
H = sparse(iH,jH,sH);
% H contains the cofficients Hei matrix
% H matrix: (nelx*nely)*(nelx*nely), sparse matrix
% The 
Hs = sum(H,2);
% Hs contain the normalization constants sigma(Hei)(i belongs to Ne)

% sum function
% S = sum(A,dim) returns the sum along dimension dim. For example, if A is 
% a matrix, then sum(A,2) is a column vector containing the sum of each row.

%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
% set the initial x equal to the prescribed volume fraction f

xPhys = x;
% Physical densities are equal to x

% The corresponding physical densities ? xe are identical to
% the design variables xe: in the sensitivity filtering approach,
% this equality always holds, while in the density filtering
% approach, it holds as long as the design variables represent
% a homogeneous field. 

% Homogeneous field: https://en.wikipedia.org/wiki/Translational_symmetry
% Compute the initial physical densities xe using equation 9 in order to
% meet the volume constraint

% For other types of filters (especially non-volume-preserving filters), it 
% may be necessary to compute the initial physical densities xe by explicit 
% application of the filter to the initial design variables xe, and to adjust
% the initial design variables in such a way that the volume constraint is 
% satisfied (as this constraint is specified in terms of the physical densities xe).

loop = 0;
change = 1;

%% START ITERATION
while change > 0.01
  % L¡Þ norm of the difference between two consecutive designs is less than 1%.
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  % sK contains the entries of the sparse stiffness matrix
  % sK depends on the physical densities x.
  % sK is obtained by reshaping the element stiffness matrix KE to obtain a
  % column vector, multiplying this vector with the appropriate Young's 
  % modulus Ee(xe) for each element and then concatenating the results for
  % all elements
  % sK = Ee(xe) * k0 = (Emin + xe^p*(E0-Emin))* k0,   k0 = KE in the
  % program
  
  K = sparse(iK,jK,sK); 
  % The assebly of global stiffness matrix K is performed by means of
  % sparse function, instead of the for-loop process, which is using in the 99 line
  % program i.e. K(edof,edof) = K(edof,edof)+x(ely,elx)^penal*KE;
  
  % The three input vectors arae iK,jK and sK.
  % sparse function: create sparse matrix
  % Syntax: S = sparse(i,j,v)
  % S = sparse(i,j,v) generates a sparse matrix S from the triplets i, j, 
  % and v such that S(i(k),j(k)) = v(k). The max(i)-by-max(j) output matrix 
  % has space allotted for length(v) nonzero elements. sparse adds together 
  % elements in v that have duplicate subscripts in i and j.
  % If the inputs i, j, and v are vectors or matrices, they must have the 
  % same number of elements. Alternatively, the argument v and/or one of the 
  % arguments i or j can be scalars.
  
  % Use repeated subscripts to accumulate values into a single sparse matrix 
  % that would otherwise require one or more loops.
  
  K = (K+K')/2;
  % The global stiffness matrix should be symmetric, because the third
  % Newton law
  % Ensure the stiffness matrix is perfectly symmetric
  % Manually cancel the rounding errors in the assembly procedure
  
  % If the stiffness matrix is sparese, symmetric and has real postive
  % diagonal elements, Cholesky factorization is used rather than LU
  % factorization, which can speed up the solving time
  
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx); % element compliance
  % efficiency improvement: compute the compliance for all elements simultaneously
  % by means of edofMat matrix
  % edofMat is used as an index into the displacement vector U to get a
  % matrix with the size of edofMat that contains the displacements
  % corresponding to the DOFs listed in the edoMat
  
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)); % sigma compliance
  
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce; % wrt the physical densities
  dv = ones(nely,nelx);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  % The application of  a sensitivity filter can be implemented as a matrix
  % product of a coefficient matrix and a vector containing the original
  % sensitivities
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:)); % r = 1e-3 
    % The sensitivities dv are identical for all elements and can therefore
    % be omitted from the definition of the heuristic updating factor
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs); % the modification of the sensitivities 
    dv(:) = H*(dv(:)./Hs); % 
  end
  
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.2;
  
  while (l2-l1)/(l1+l2) > 1e-3 % Stop condition is specified in relative terms 
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));  
    % sqrt(-dc./dv/lmid) equals Be
    % Bisection to search ¦Ë, 
    % the sensitivity dv of the volume constraint is explicitly taken into account
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs; % density filtering
    end
    if sum(xPhys(:)) > volfrac*nelx*nely %Lagrange multiplier lmid is 
        % Determine the Lagrange multiplier using the physical densities 
        % instead of the design variables
        
        % The density filter is volume-preserving so that the volume
        % constraint can equally well be evaluated in terms of the design
        % variables
        
        % If the filter is non-volume preserving, it is necessary to
        % evaluate the volume constraint in terms of the physical densities
        l1 = lmid;
    else
        l2 = lmid;
    end
  end
  % The termination means the Lagrangian multiplier is converge
  change = max(abs(xnew(:)-x(:)));
  % change: L¡Þ norm of the difference in terms of design variables
  x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

