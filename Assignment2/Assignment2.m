%% MESF6910 Topology Optimization  
%% Assignment #2 FEM Analysis of a Michell-Type Structure---Stess_Strain_Analysis  
%% Source Code

%%% Revised by CHEN KAIWEN 
%%% Modify the Finite_Element Code,state variables and boundary conditions
%%% Based on "A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLESIGMUND, OCTOBER 1999"
%%% Problem Description: Assignment 2 Strain-Stress.ppt
%% Assumption

% Design domain is assumed to be rectangular 12m * 6m
% Discretized by quadrilateral finite elements (0.1m * 0.1m) 
% Using Bilinear Rectangle plane element: The distribution of displacement 
% within the element is linear along x and y axis
% The displacement u at any interior point of an element is interpolated in
% terms of the shape functions N and the nodal displacements d as
% u=N*d
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
% Poisson＊s ratio = 0.3
%% Main Program
% nelx=120
% nely=60
% volfrac=1
% penal=1
% Both nodes and elements are numbered column wise from left to right
% the numbering of elements and nodes is column by column starting in the 
% upper left corner)
% ----------------------------------------------------------------------


[VMStress,VMStrain,Straintotal,Stresstotal] = Sress_Stain(120,60,1);
f1=figure('Name','VonMisesStress','NumberTitle','off');
f1=imagesc(VMStress);
colorbar;

f2=figure('Name','VonMisesStrain','NumberTitle','off');
f2=imagesc(VMStrain);
colorbar;

function [VonMisesStress,VonMisesStrain,Strain,Stress] = Sress_Stain(nelx,nely,penal)
syms b h y x
% INITIALIZE
Strainleftbottom = cell(nely,nelx);
Strainleftupper = cell(nely,nelx);
Strainrightupper = cell(nely,nelx);
Strainrightbottom = cell(nely,nelx);
Strain = cell(nely,nelx,2,2);
Straincenter = cell(nely,nelx);


Stressleftbottom = cell(nely,nelx);
Stressleftupper = cell(nely,nelx);
Stressrightupper = cell(nely,nelx);
Stressrightbottom = cell(nely,nelx);
Stress = cell(nely,nelx,2,2);
Stresscenter = cell (nely,nelx); 
Deviatoriastrain = cell (nely,nelx);
VonMisesStress = zeros(nely,nelx);
VonMisesStrain = zeros(nely,nelx);

B = 1/(4*b*h) * [-(h-y) 0 (h-y) 0 (h+y) 0 -(h+y) 0; % B matrix: The relationship between strain and displacement
               0 -(b-x) 0 -(b+x) 0 (b+x) 0 (b-x);
               -(b-x) -(h-y) -(b+x) (h-y) (b+x) (h+y) (b-x) -(h+y)];
Bleftbottom = eval(subs(B,[b,h,y,x],[0.05,0.05,-0.05,-0.05]));  % B matrix at the left bottom corner
Bleftupper = eval(subs(B,[b,h,y,x],[0.05,0.05,0.05,-0.05])); % b= 0.05m h=0.05m 
Brightupper = eval(subs(B,[b,h,y,x],[0.05,0.05,0.05,0.05])); 
Brightbottom = eval(subs(B,[b,h,y,x],[0.05,0.05,-0.05,0.05])); 
Bcenter = eval(subs(B,[b,h,y,x],[0.05,0.05,0,0]));

% Material Constants 
% E means Young's modulus(modulus of elasticity)  100 Mpa
E = 100000000; % unit:Pa
nu = 0.3; % nu means Poisson＊s ratio
 
% Plane Stress Stiffness matrix
StiffnessMatrix = E/(1-nu^2).*[1 nu 0
                              nu 1 0
                              0 0 (1-nu)/2];

% FE-ANALYSIS
% Displacement calculation
[U]=FE(nelx,nely,penal);
Displacement=U(:,1)+U(:,2)+U(:,3); % Three forces case

% Stress-Strain Analysis
% Note: The node displacement order is given by the method using in the Assignment 2 Strain-Stress.ppt
% Anti-clockwise, not the same as in the displacement analysis.
for ely = 1:nely
    for elx = 1:nelx
        n1 = (nely+1)*(elx-1)+ely; % variable n1:upper left element node numbers in global node numbers
        n2 = (nely+1)* elx +ely; % variable n2:upper right elemnt node number in global node numbers
        u1 = Displacement(2*n1+1);
        v1 = Displacement(2*n1+2);
        u2 = Displacement(2*n2+1);
        v2 = Displacement(2*n2+2);
        u3 = Displacement(2*n2-1);
        v3 = Displacement(2*n2);
        u4 = Displacement(2*n1-1);
        v4 = Displacement(2*n1);
        nodaldisplacement=[u1;v1;u2;v2;u3;v3;u4;v4];
        
        % Strain Matrix=[exx;eyy;2exy] 
        % In the Strain Matrix, the third element is 2exy=rxy, so it should
        % be {exx;eyy;rxy}
        Strainleftbottom{ely,elx} = Bleftbottom * nodaldisplacement;
        Strainleftupper{ely,elx} = Bleftupper * nodaldisplacement;
        Strainrightupper{ely,elx}= Brightupper * nodaldisplacement;
        Strainrightbottom{ely,elx}= Brightbottom * nodaldisplacement;
        Straincenter{ely,elx} = Bcenter * nodaldisplacement;
        
        Strain{ely,elx,1,1} = Strainleftupper{ely,elx}; 
        Strain{ely,elx,1,2} = Strainrightupper{ely,elx}; 
        Strain{ely,elx,2,1} = Strainleftbottom{ely,elx};
        Strain{ely,elx,2,2} = Strainrightbottom{ely,elx};
        
        % https://dianafea.com/manuals/d944/Analys/node405.html
        Deviatoriastrain{ely,elx} = [ 2/3 -1/3 -1/3; -1/3 2/3 -1/3; -1/3 -1/3 2/3]...
                        * [Straincenter{ely,elx}(1);Straincenter{ely,elx}(2);0];
                    
        VonMisesStrain(ely,elx) = 2/3 * sqrt(3 * (Deviatoriastrain{ely,elx}(1)^2 ...
        + Deviatoriastrain{ely,elx}(2)^2 + Deviatoriastrain{ely,elx}(3)^2) / 2 + ...
        3 * (Straincenter{ely,elx}(3)^2)/4);
                                 
        
        % Hooke's Law for Plane Stress 
        % http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_stress.cfm
        % Here is using the Plane Stress Hooke's Law via Engineering Shear Strain
        % rxy is the same as 2exy. G means shear modulus. G = E/(2(1+v))
        % 考xy= G * rxy(2exy)
        % StressMatrix = {考xx;考yy;考xy}
        Stressleftbottom{ely,elx} = StiffnessMatrix * Strainleftbottom{ely,elx};
        Stressleftupper{ely,elx} = StiffnessMatrix * Strainleftupper{ely,elx};
        Stressrightupper{ely,elx} = StiffnessMatrix * Strainrightupper{ely,elx};
        Stressrightbottom{ely,elx} = StiffnessMatrix * Strainrightbottom{ely,elx};
        Stress{ely,elx,1,1} = Stressleftupper{ely,elx}; 
        Stress{ely,elx,1,2} = Stressrightupper{ely,elx}; 
        Stress{ely,elx,2,1} = Stressleftbottom{ely,elx};
        Stress{ely,elx,2,2} = Stressrightbottom{ely,elx};
        
        Stresscenter{ely,elx} = StiffnessMatrix * Straincenter{ely,elx};
        
        
        % General Plane Stress
        % https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
        VonMisesStress(ely,elx) = sqrt(Stresscenter{ely,elx}(1)^2 ...
            - Stresscenter{ely,elx}(1) * Stresscenter{ely,elx}(2) + ...
            Stresscenter{ely,elx}(2)^2 + 3 * Stresscenter{ely,elx}(3)^2);

    end
end



end




%% Finite_Element Code
function [U]=FE(nelx,nely,penal)
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
        K(edof,edof) = K(edof,edof)+1^penal*KE; % 1=x which equals the volumefraction 
    end
 end

% DEFINE LOADS AND SUPPORTS
 F(3782,1) =-15;% Apply P2 at the right buttom corner of element(elx=30 ely=60) 
 % 2*(nely+1)=122˙122 is the y direction of left buttom corner node in
 % global
 % 122*31=3782˙122*61=7442˙122*91=11102;
 F(7442,2) =-30;% Apply P1 at the buttom right corner of element(elx=60 ely=60)
 F(11102,3)=-15;% Apply P2 at the right buttom corner of element(elx=90 ely=60)
 fixeddofs =[121 122 14762];% fixeddofs define the degrees of freedom that are fixed
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
 E = 100000000; % unit:Pa
 % nu means Poisson＊s ratio
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
