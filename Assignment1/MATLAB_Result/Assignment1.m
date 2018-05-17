%% MESF6910 Topology Optimization 
%% Assignment #1 FEM Analysis of a Michell-Type Structure
%% Source Code

%%% Revised by CHEN KAIWEN 
%%% Modify the Finite_Element Code,state variables and boundary conditions
%%% Based on "A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLESIGMUND, OCTOBER 1999"
%%% Problem Description: Assignment 1.docx 
%% Assumption

% Design domain is assumed to be rectangular 12m * 6m
% Discretized by quadrilateral finite elements (0.1m * 0.1m) 
%% Glossary

% volfrac: Volume fraction
% penal: pnealization power
% nelx: number of elements in the horizontal direction 
% nely: number of elements in vertical direction
% nelx*nely: number of elements used to discretize the design domain
% rmin: round,filter size(divided by element size),not used in FE
% x: the vector of design variables 
% xmin: a vector of minimum relative densities
%% Settings 

% set volfrac and penal equal to 1
% set nelx equals to 120
% set nely equals to 60
% modulus of elasticity = 100MPa
% Poisson¡¯s ratio = 0.3
%% Main Program
% nelx=120
% nely=60
% volfrac=1
% penal=1
% Both nodes and elements are numbered column wise from left to right
% the numbering of elements and nodes is column by column starting in the 
% upper left corner)
% ----------------------------------------------------------------------
function top(nelx,nely,volfrac,penal)

% INITIALIZE
x(1:nely,1:nelx) = volfrac;

% FE-ANALYSIS
[U]=FE(nelx,nely,x,penal);
Displacement=U(:,1)+U(:,2)+U(:,3); % Three forces case

end




%% Finite_Element Code
function [U]=FE(nelx,nely,x,penal)
 [KE] = lk;
 K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1)); % K is global stiffness matrix
 F = sparse(2*(nely+1)*(nelx+1),3);% F is global force vector   
 U =sparse(2*(nely+1)*(nelx+1),3); % U is global displacement vector  
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

% DEFINE LOADS AND SUPPORTS
 F(3782,1) =-15;% Apply P2 at the right buttom corner of element(elx=30 ely=60) 
 % 2*(nely+1)=122£»122 is the y direction of left buttom corner node in
 % global
 % 122*31=3782£»122*61=7442£»122*91=11102;
 F(7442,2) =-30;% Apply P1 at the buttom right corner of element(elx=60 ely=60)
 F(11102,3)=-15;% Apply P2 at the right buttom corner of element(elx=90 ely=60)
 fixeddofs =[121 122 14762];% fixeddofs define the degrees of freedom that are fixed
 alldofs = [1:2*(nely+1)*(nelx+1)];
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
 E = 100000000; % unit:Pa
 % nu means Poisson¡¯s ratio
 nu = 0.3;
 % Integral of element level
 % Analytical solution for k and KE
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
 end
end