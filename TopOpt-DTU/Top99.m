%% MESF6910 Topology Optimization_A 99 LINE TOPOLOGY OPTIMIZATION CODE
%% Source Code

%%% Based on "A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLESIGMUND, OCTOBER 1999"
%% Assumption

% Design domain is assumed to be rectangular 
% Discretized by square finite elements 
%% Glossary

% volfrac: Volume fraction
% penal: pnealization power
% nelx: ratio of elements in the horizontal direction 
% nely: ratio of elements in vertical direction
% nelx*nely: number of elements used to discretize the design domain
% rmin: round,filter size(divided by element size)
% x: the vector of design variables
% xmin: a vector of minimum relative densities
%% 
%% Main Program
% ----------------------------------------------------------------------
function Top99(nelx,nely,volfrac,penal,rmin)

% INITIALIZE
x(1:nely,1:nelx) = volfrac;
loop = 0;
change = 1.;

% START ITERATION
while change > 0.01  % change = max(max(abs(x-xold))) in loop;
loop = loop + 1;
xold = x;

% FE-ANALYSIS
[U]=FE(nelx,nely,x,penal);

% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
[KE] = lk; % [KE] = k0 in the paper
c = 0.;
for ely = 1:nely
 for elx = 1:nelx
  n1 = (nely+1)*(elx-1)+ely;
  n2 = (nely+1)* elx +ely;
  Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
  c = c + x(ely,elx)^penal*Ue'*KE*Ue; % c:compliance, all all elements together
  dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue; % partial c partial xe
 end
end

% FILTERING OF SENSITIVITIES
[dc] = check(nelx,nely,rmin,x,dc); % equation (5) in the Paper

% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
[x] = OC(nelx,nely,x,volfrac,dc);

% PRINT RESULTS
change = max(max(abs(x-xold)));
disp(['It.:' sprintf('%4i',loop) 'Obj.:' sprintf('%10.4f',c) ...
' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
' ch.: ' sprintf('%6.3f',change )])

% PLOT DENSITIES
colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
end

%% Optimality Criteria based optimizer 
function [xnew]=OC(nelx,nely,x,volfrac,dc)
 l1 = 0; l2 = 100000; move = 0.2;
 while (l2-l1 > 1e-4)
 lmid = 0.5*(l2+l1);
 xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
if sum(sum(xnew)) - volfrac*nelx*nely > 0
    l1 = lmid;
else
 l2 = lmid;
end
 end

%% MESH-INDEPENDENCY FILTER
% Matlab implementation of modifying the element sensitivities,equation (5)
% in the paper

% Set the rmin greater than or equal to 1 so that use the filter! 
% If rmin set less than 1, the filtered sensitivities will be equal to the
% original sensitivities making the filter inactive.

% Only the elements within a saqure with side lengths 
% 2 * round(rmin) around the considered element are searched in the design domain 
% in order to find the elements that lie within the radius rmin

function [dcn]=check(nelx,nely,rmin,x,dc)
 dcn=zeros(nely,nelx);
 for i = 1:nelx
     for j = 1:nely
        sum=0.0;
            for k = max(i-round(rmin),1):min(i+round(rmin),nelx)
                for l = max(j-round(rmin),1):min(j+round(rmin), nely)
                    fac = rmin-sqrt((i-k)^2+(j-l)^2); % fac = Hf = rmin - distance(e,f)
                    sum = sum+max(0,fac); % Hf is zero outside the filter area , if Hf is smaller than 0, use 0 instead
                    dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k); % x(l,k) = xf in Paper, dc(l,k) = partial c partial xf
                end
            end
        dcn(j,i) = dcn(j,i)/(x(j,i)*sum); % sum = sigma Hf, dcn = <partial c partial xe> Hat 
     end
 end

%% Finite_Element Code
function [U]=FE(nelx,nely,x,penal)
 [KE] = lk;
 K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1)); % K is global stiffness matrix
 F = sparse(2*(nely+1)*(nelx+1),1);% F is global force vector   
 U =sparse(2*(nely+1)*(nelx+1),1); % U is global displacement vector  
 for ely = 1:nely
    for elx = 1:nelx
        % insert the element stiffness matrix at the right place in the
        % global stiffness matrix
        n1 = (nely+1)*(elx-1)+ely; % variable n1:upper left element node numbers in global node numbers
        n2 = (nely+1)* elx +ely; % variable n2:upper right elemnt node number in global node numbers
        edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1;2*n2+2;2*n1+1; 2*n1+2];
        K(edof,edof) = K(edof,edof)+x(ely,elx)^penal*KE;
    end
 end

% DEFINE LOADS AND SUPPORTS(HALF MBB-BEAM)
 F(2,1) = -1; % Apply a vertical unit force in the upper left corner       
 fixeddofs = union(1:2:2*(nely+1),2*(nelx+1)*(nely+1));% fixeddofs define the degrees of freedom that are fixed
 alldofs = 1:2*(nely+1)*(nelx+1);
 freedofs = setdiff(alldofs,fixeddofs);% freedofs indicate the degrees of freedom which are unconstrained
 
 % setdiff operator finds the free degrees of freedoms as the difference 
 % between all degrees of freedom and the fixed degrees of freedom
 
 
 % SOLVING
 % Implement supports, by eliminating fixed degrees of freedom from the
 % linear equations
 U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:); 
 U(fixeddofs,:)= 0;
 
 %%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%
 % calculate element stiffness matrix 
 % 8 by 8 matrix for a square bi-linear 4 node element
function [KE]=lk
 % E means Young's modulus(modulus of elasticity)  100 Mpa
 E = 1;
 % nu means Poisson¡¯s ratio
 nu = 0.3;
 k=[ 1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ...
 -1/4+nu/12 -1/8-nu/8 nu/6 1/8-3*nu/8];
 KE = E/(1-nu^2)* [ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
 k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
 k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
 k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
 k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
 k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
 k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
 k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];