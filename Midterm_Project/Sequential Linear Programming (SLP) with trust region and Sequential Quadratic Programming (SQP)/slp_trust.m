function [x,f,Converged,output,lambda]=slp_trust(Fun,X0,Options,vlb,vub,Grd,varargin)
% Sequential Linear Programming (SLP) with a Trust Region Strategy
% finds the constrained minimum of a function of several variables.
% Version 1.3
%
% Copyright (c) 2015, Robert A. Canfield. All rights reserved.
%                     See accompanying LICENSE.txt file for conditions.
%
%  fmincon-compatible problem structure input argument
%          optimtool GUI option "Export to Workspace" dialog box 
%          sends problem information to the MATLAB workspace as a structure
%
%          usage: [x,f,Converged,output,lambda]=slp_trust( problem )
%
%          input: problem   - Data structure with fields:
%                 objective - Objective function
%                 x0        - Initial point for x
%                 Aineq     - Matrix for linear inequality constraints
%                 bineq     - Vector for linear inequality constraints
%                 Aeq       - Matrix for linear equality constraints
%                 beq       - Vector for linear equality constraints
%                 lb	       - Vector of lower bounds
%                 ub        - Vector of upper bounds
%                 nonlcon   - Nonlinear constraint function 
%                 options   - Options created with optimset
%
%  Optimization Toolbox Version 1-compatible input arguments
%
%  usage: [x,f,Converged,output,lambda]=slp_trust(Fun,X0,Opts,vlb,vub,Grd,P1,P2,...)
%
%  input:   Fun     - function handle (or string) which returns the
%                     value of the objective function and a vector of
%                     constraints, i.e., [f,g]=fun(x,active,P1,P2,...) 
%                     where x is the design variable vector and active is
%                     an index/boolean vector pointing to active constraints.
%                     f is minimized such that g<=zeros(g).
%           X0      - initial vector of design variables
%           Options - (optional) a structure according to optimset & fields
%                     TrustRegion: on, off, simple, merit, TRAM (string)
%                     ComplexStep: off, on (string)
%                     cutg:        active constraint cutoff value (real)
%           vlb     - (optional) vector of lower bounds on design variables
%           vub     - (optional) vector of upper bounds on design variables
%           Grd     - (optional) function handle that returns a vector of 
%                     function gradients and a matrix of constraint gradients
%                     i.e. [fp,gp]=grd(x,active,P1,P2,...) where
%                     active is an index (or boolean) vector to active g.
%           Pn      - (optional) variables directly passed to Fun and Grd
%
%           Note: optional inputs can be skipped by inputing []
%
%  output: x         - vector of design variables at the optimal solution
%          f         - final value of objective function
%          Converged - Convegence flag
%          output   - Structure of output results (iteration history)
%          lambda   - Structure of Lagrange multipliers
%
%  Written by:    Robert A. Canfield
%  e-mail:        bob.canfield@vt.edu
%
%  Created:        4/13/06
%  Last modified:  4/18/17
%
% The function format is based on the MATLAB function constr.m written
% by Andy Grace of MathWorks, 7/90, supplemented with optimization toolbox 2
% structures for options and output.
%
% Trust Region (TR) algorithm and filter based on: 
% Nocedal and Wright, Numerical Optimization. New York: Springer, 2006 
% Algorithms 4.1 for TR & 15.1 for filter; Equation (15.4) for simple merit
%--------------------------------------------------------------------------
% Copyright (c) 2015, Robert A. Canfield. All rights reserved.
%                     Virginia Tech
%                     bob.canfield@vt.edu
%                    <http://www.aoe.vt.edu/people/faculty/canfield.html>
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal with the Software without restriction, including without 
% limitation the rights to use, copy, modify, merge, publish, distribute, 
% sublicense, and/or sell copies of the Software, and to permit persons 
% to whom the Software is furnished to do so, subject to the following 
% conditions:
% 
% * Redistributions of source code must retain the above copyright notice,
%   this list of conditions and the following disclaimers.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimers in the
%   documentation and/or other materials provided with the distribution.
% 
% * Neither the names of Robert A. Canfield, Virginia Tech, 
%   Air Force Institute of Technology, nor the names of its contributors 
%   may be used to endorse or promote products derived from this Software 
%   without specific prior written permission.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
% OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
% THE USE OR OTHER DEALINGS WITH THE SOFTWARE. 
%--------------------------------------------------------------------------

%--Modifications
%  10/30/11 - fdgrd internal function added
%   9/27/15 - active constraint strategy added
%   9/29/15 - fmincon compatibility with single input: problem structure
%   9/30/15 - fix grd -> eval grdstr
%  10/16/15 - user Fcn, Grd with or without active
%  11/12/15 - exit if checkbounds msg error, reset X0 if needed
%  12/09/15 - FinDiffRelStep not backward compatible
%   9/28/16 - gradf=gradg typo
%   10/2/16 - refine bindLM
%   11/7/16 - TR.Bound returned
%   11/8/16 - active0
%  11/12/16 - display trust region filter status string; lp optimoptions
%  11/13/16 - return all output arguments from trust_region in TR
%  11/14/16 - Active constraint set for gradients; TR.cutg
%  11/24/16 - linprog no feasible solution, find least infeasible
%  11/25/16 - R2016b optim toolbox version 7.5 optimoptions compatibility
%   3/21/17 - Detect user unaware of active when varargin supplied
%   4/01/17 - MPEA
%   4/ 9/17 - nac<ncon fixed
%   4/18/17 - Use Merit in Slowed only for Trust

%% Check inputs
if nargin<1
   disp('usage: [x,f,exitflag,output,lambda]=slp_trust(fun,x0,options,vlb,vub,grd,P1,P2,...)')
   return
end
if nargin>1
   if nargin<3, Options=[];end
   if nargin<4, vlb=[]; end
   if nargin<5, vub=[]; end
   if nargin<6, Grd=[]; end; fd_gradients = isempty(Grd);
else
   Problem = Fun;
   [Obj,X0,A,b,Aeq,beq,vlb,vub,Con,Options]=sqpcrkfmincon(Problem);
   Fun = @(x) fun2(x,Obj,Con,A,b,Aeq,beq);
   Grd = @(x) grd2(x,Obj,Con,A,Aeq);
   fd_gradients = nargout(Obj)<2 || nargout(Con)<3;
end

%% Process options
% Standard optim toolbox options
if  isempty(Options)
   options = optimset('fmincon','Display','iter');
   if isempty(options.TolX), options.TolX=options.TolCon; end
   TolOpt  = options.TolFun / max(1,optionsTolX);
elseif isstruct(Options)
   TR = Options;
   options = optimset(optimset('fmincon'),Options);
   if isempty(options.TolX), options.TolX = options.TolCon; end
   if isfield(Options,'OptimalityTolerance')
      TolOpt = Options.OptimalityTolerance;
   elseif isfield(Options,'TolOpt')
      TolOpt = Options.TolOpt;
   else
      TolOpt = options.TolFun / min(1,options.TolX);
   end
elseif isa(Options,'optim.options.Fmincon')
   warning('off','MATLAB:structOnObject');
   options = struct(Options);
   TolOpt = options.OptimalityTolerance;
   options.TolFun = TolOpt;
   options.TolCon = options.ConstraintTolerance;
   options.TolX   = options.ConstraintTolerance;
else
   error('slp_trust:options','Options is not viable structure')
end
TolFun = options.TolFun;
TolCon = options.TolCon;
TolX   = options.TolX;
Tolvb  = optimget(optimset('fmincon'),'TolCon'); % tolerance for active DV
if isempty(options.MaxIter), options.MaxIter = 50*length(X0);  end
% Supplemental options for slp_trust
if isfield(options,'ComplexStep') && strcmpi(options.ComplexStep,'on')
   options.FinDiffType = 'complex'; 
end
if ~isfield(options,'cutg')   || isempty(options.cutg),   options.cutg=0.1; end
if ~isfield(options,'manyDV') || isempty(options.manyDV), options.manyDV=50; end
if ~isfield(options,'manyCon')|| isempty(options.manyCon),options.manyCon=50;end
TR.options = options;
MPEA = isfield(TR,'MPEA') && strcmpi(TR.MPEA,'on');

%% Check lower and upper bounds: vlb and vub
ndv = length(X0(:)); % number of design variables
[x,vlb,vub,msg] = checkbounds(X0,vlb,vub,ndv); % Optimization toolbox routine
if ~isempty(msg), error(msg), end
lenvlb=length(vlb); ilb=(1:lenvlb)';
lenvub=length(vub); iub=(1:lenvub)';
if lenvlb && any(x(ilb)<vlb)
   X0(ilb)=max(x(ilb),vlb);
   warning('Initial x vector was not within lower bounds.');
   disp('Reset x to '), disp(X0)
end
if lenvub && any(x(iub)>vub)
   X0(iub)=min(x(iub),vub);
   warning('Initial x vector was not within upper bounds.');
   disp('Reset x to '), disp(X0)
end

%% Set up function and gradient calls. Initial evaluations.
if all(size(x)==size(X0))
   xshape = @(x) x;
else
   xshape = @(x) reshape(x(:),size(X0,1),size(X0,2));
end;
if nargin(Fun) > 1+numel(varargin)
   fun   = @(x)    FunEval(xshape(x));
   funfd = @(x,active) Fun(xshape(x),active,varargin{:});
   active=NaN; % flag to evaluate all constraints (first time)
elseif nargin(Fun) == 1+numel(varargin)
   fun   = @(x)        Fun(xshape(x),varargin{:});
   funfd = @(x,active) Fun(xshape(x),varargin{:}); 
else
   fun   = @(x) Fun(xshape(x));
   funfd = fun; 
end
if fd_gradients
   grdstr = 'fdgrd(funfd,x,f,g,options,active)';
elseif nargin(Grd) > 1+numel(varargin)
   grdstr = 'Grd(xshape(x),active,varargin{:})';
elseif nargin(Grd) == 1+numel(varargin)
   grdstr = 'Grd(xshape(x),varargin{:})';
else
   grdstr = 'Grd(xshape(x))';
end

% Initial function evaluations.
[f, g]  = fun(x); g=g(:);
[mg,mj] = max(g);
ncon    = length(g);
if nargin(Fun)<2+numel(varargin) % user unaware of active constraints
   options.cutg = Inf;
   TR.cutoff    = Inf;
   nretain = ncon;
else
   nretain = 2*ndv; % Retain more active constraints than design variables
end
[gradf,gradg,active,cutoff] = GradEval(x,[],nretain);

if length(gradf) ~= ndv
   error('slp_trust:gradf','Objective gradient length does not match # design variables')
elseif size(gradg,1) ~= ndv
   error('slp_trust:gradg','# rows of gradg should be # design variables')
% elseif size(gradg,2) ~= ncon % but maybe nretain,nac<ncon
%    error('slp_trust:nac','# columns of gradg should be # constraints')
elseif nargin(Fun)>1+numel(varargin) % user aware of active constraints
   TR.cutoff = cutoff;
   nretain = ndv; % Subsequently retain at least ndv active constraints
end
if nargout>3
   output.f(1) = f;
   output.g(:,1) = g;
end

% Set Trust Region (TR) data structure fields
% TR.SimpleTrust=~isfield(TR,'TrustRegion') || ... % turn off SimpleTrust
%                 isempty(TR.TrustRegion)   || ... % as default
%                ~isempty(strfind(lower(TR.TrustRegion),'simple'));
if MPEA
   TR.gradf=gradf; TR.gradg=gradg; TR.f=f; TR.g=g; TR.x=x; 
   TR.p=ones(size(x)); TR.iter=0;
end
TR = trust_region(TR,x,[],vlb,vub,f,f,gradf,g(active),g(active),gradg);

%% SLP loop
% ------------------------------------------------------------------------
Iter=0;
TR.stat=0;
Converged=false;
Slowed=false;
if ndv>options.manyDV
   optLP = optimoptions(@linprog,'Display','off','Algorithm','interior-point');
elseif verLessThan('optim','7.5')
   optLP = optimoptions(@linprog,'Display','off','Algorithm','active-set'); %#ok<LINPROG>
   warning('off','optim:linprog:AlgOptsWillError');
elseif mg<TolCon && ... % dual-simplex won't return best feasible solution
       any(gradg(:)~=0) % crashes when no constraints (g<0,gradg=0)
   optLP = optimoptions(@linprog,'Display','off','Algorithm','dual-simplex');
else % interior-point-legacy will return least infeasible design
   optLP = optimoptions(@linprog,'Display','off','Algorithm','interior-point-legacy');
end
if strcmpi(options.Display,'iter')
   indent=blanks(9); if MPEA, indent=[indent,'MPEA ']; end
   disp(' ')
   disp([indent, 'Sequential Linear Programming Iteration History'])
   disp('Iteration      Objective MaxConstraint    Index   Step-size   Merit      MoveLimit  TrustRatio')
   fprintf('%9.0f %14.5g  %12.4g  %5.0f  %10.4g  %10.4g\n',[Iter f mg mj 0 TR.Merit])
end
%                                        (allow LP infeasible to continue)
while ~Converged && Iter<options.MaxIter && (TR.stat>=0 || TR.stat<=-2)
   Iter = Iter + 1;
   x0 = x;
   f0 = f;
   g0 = g;
   mg0 = mg;
   active0 = active;
   MoveLimit0 = TR.MoveLimit;
   dxlb = max( vlb-x0, -TR.delx );
   dxub = min( vub-x0,  TR.delx );
   b = g0(active);
   if MPEA && Iter > 1
      [dx,TR,out,lambda]=mpeaprog(gradf,gradg,f,b,lambda,x,dxlb,dxub,TR,optLP);
   else
      [dx,~,TR.stat,out,lambda]=linprog(gradf,gradg',-b,[],[],dxlb,dxub,[],optLP);
      if isempty(dx) && (TR.stat==-2 || TR.stat==-5) % No feasible point found.
         j  = g(active)>TolCon;                              % Add slack.
         c  = [gradf+gradg(:,j)*max(TR.penalty(j),1);1]; % merit fcn is obj.
         gl = g(active).*max(TR.penalty,1); % weighted linearized constraints
         A  = [(gradg*diag(max(TR.penalty,1))).', -gl];
         [dx,~,TR.stat,out,lambda]=linprog(c,A,-gl,[],[],[dxlb;0],[dxub;1],[],optLP);
         dx   = dx(1:end-1);
         gapx = g0(active) + gradg'*dx;
         j    = gapx>TolCon;
         lambda.ineqlin(j) = lambda.ineqlin(j) + gl(j)*lambda.upper(end);
         lambda.lower      = lambda.lower(1:end-1);
         lambda.upper      = lambda.upper(1:end-1);
      end
      TR.fapx = f0         + gradf'*dx;
      TR.gapx = g0(active) + gradg'*dx;
   end
   x = max(vlb, min(vub, x0 + dx));
   [f,g] = fun(x); g=g(:);
   [mg,mj] = max(g);

%% Update trust region move limits.
   TR = trust_region(TR,x,dx,vlb,vub,f0,f,gradf,g0(active),g(active),gradg,lambda);
   active1 = TR.bindLM;

   % Check convergence.
   if ~TR.rejected
      [gradf,gradg] = GradEval(x,active,nretain); % gradients
      Feasible = mg <= TolCon;
      Slowed   = abs(f-f0) < TolFun && max(abs(dx)) < TolX;
      if TR.trust && ~TR.filter.flag
         Slowed = Slowed && abs(f-TR.Merit) < TolFun;
      end
      bound  = x<vlb+Tolvb | x>vub-Tolvb; % active side constraints
      Lagrangian  = gradf(:) + (lambda.upper - lambda.lower).*(bound);
      if ~isempty(TR.bindLM) % recompute Lagrange multpliers
         if Slowed && any(~bound) && any(lambda.ineqlin(TR.bindLM))
            Lambda = -gradg(~bound,TR.bindLM) \ gradf(~bound);
            if any(Lambda<0), Lambda=lambda.ineqlin(TR.bindLM); end
         else
            Lambda=lambda.ineqlin(TR.bindLM);
         end
         Lagrangian = Lagrangian + gradg(:,TR.bindLM)*Lambda;
      end
      KKT = Feasible && norm(Lagrangian,inf) < TolOpt ...
                     && max(abs(dx))<TolX;
      Converged = KKT || (Slowed && Feasible);
   else % if rejected
      newactive = ~ismember(active1,active0);
      if any(newactive)
         active = active1(newactive);
         [~,gradg] = GradEval(x,active);
         active = union(active0,active1); %#ok<*NASGU>
      end
      Lagrangian = NaN;
   end

   % Output
   if strcmpi(options.Display,'iter') % Print intermediate results
      if Converged || Slowed
         if TR.Bound, BoundStr=' Bound'; else, BoundStr=' Unbound'; end
         TRstatus = strcat(TR.filter.str,BoundStr);
      else
         TRstatus = TR.filter.str;
      end
      disp([sprintf('%9.0f %14.5g  %12.4g  %5.0f  %10.4g  %10.4g  %9.4g  %10.4g  ',...
           [Iter f mg mj max(abs(dx)) TR.Merit MoveLimit0 TR.Ratio]) TRstatus])
   end
   if nargout>3
      output.f(end+1) = f;
      output.g(:,end+1) = g;
   end
   if TR.rejected
      x=x0; f=f0; g=g0; mg=mg0; % Reset to previous point
   end
end
if nargout>3, output.iterations = Iter; output.message=TR.stat; end
if nargout>4, lambda.ineq = lambda.ineqlin; end
%
% ------------------------------------------------------------------------

%% Print final results
if strcmpi(options.Display,'iter')
   disp('              ----------  ------------         ----------')
   fprintf('    Criteria   %9.4g  %12.4g         %10.4g\n',[TolFun TolCon TolX])
end
if ~strcmpi(options.Display,'off')
   if Converged
      if KKT
         stopped='converged. '; 
      elseif Slowed && Feasible
         stopped='slowed.    ';
      else
         stopped='           ';
      end
      disp(['SLP ',stopped,'Final objective function value = ',num2str(f)])
      disp(['               Lagrangian gradient   2-norm = ',num2str(norm(Lagrangian))])
      disp(['               Lagrangian gradient inf-norm = ',num2str(norm(Lagrangian,inf))])
      disp(['               Optimality Tolerance         = ',num2str(TolOpt)])
      TolLM = sqrt(eps)*max([1;lambda.ineqlin]);
      active = find(lambda.ineqlin>TolLM);
      if strcmp(options.Display,'Iter') && numel(active)<=20
         disp( '               Lagrange Multipliers   (j)')
         if any(active)
            fprintf('%34.4g  %4.0f\n',[lambda.ineqlin(active), active]');
         end
         bound = find( lambda.lower>Tolvb | lambda.upper>Tolvb );
         if any(bound)
             disp( '               Lower    Upper         (i)')
             fprintf('%22.4g  %10.4g  %4.0f\n',...
                 [lambda.lower(bound),...
                  lambda.upper(bound), bound]')
         end
      end
   else
      disp(['SLP did NOT converge in ',num2str(Iter),' iterations.'])
   end
   if TR.trust
      if TR.filter.flag
         merit_filter = 'filter';
      else
         merit_filter = 'Merit function';
      end
      if TR.SimpleTrust
         simple = 'simple ';
      else
         simple=[];
      end
      disp(['Trust Region Strategy uses ',simple,merit_filter])
   end
   disp('* Dominates prior points')
   disp('+ Nondominated')
   disp('- Dominated by prior point(s)')
   if TR.trust && TR.adapt
      disp('! Trust Radius set by Merit function minimization')
      disp('_ Trust Radius set by target Trust Ratio')
   end
   if TR.filter.flag || TR.adapt
      disp('f/g/m Objective/Constraint/Merit governs Trust Ratio')
   end
end
if TR.stat<0 && ~Converged
   error('slp_trust:linprog',out.message)
end


%% Nested internal functions
   function [f,g]=FunEval(x)
      if any(isnan(active)) % flag to evaluate all constraints first time
         [f,g] = Fun(xshape(x),active,varargin{:});
      else % evaluate all constraints each outer-loop iteration
         [f,g] = Fun(xshape(x),(1:ncon),varargin{:});
      end
      g=g(:);
   end

   function [gradf,gradg,active,cutoff] = GradEval(~,active,nretain)
      if nargin<2, active=[]; end
      if nargin<2 || (nargin>2 && nretain>=ncon) % evaluate all gradients
         cutoff = inf;
         active = 1:ncon; % active used in grdstr
      elseif isempty(active)
         cutoff = max([abs(options.cutg),max(mg),10*TolCon]);
         active  = find(g>=min(0,mg)-cutoff);
      else
         cutoff = abs(min(g(active)));
      end
      if nargin>2 && nretain<ncon
         [~,gorder]=sort(g,'descend');
         cutoff = max(cutoff,abs(g(gorder(nretain))));
         active  = find(g>=min(0,mg)-cutoff);
      end
      [gradf,gradg] = eval(grdstr);  gradf=gradf(:);
   end

end
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Internal functions
%
%% Finite Difference Gradient sub-function
function [df,dg] = fdgrd( fg, x0, f0, g0, options, active ) %#ok<DEFNU>
%  FDGRD      Calculates first forward finite difference gradients
%             of objective and constraint functions.
%
%  usage:        [df,dg]=sqpfdgrd(fcn, x0, f0, g0, xmin, xmax, x1,P1,...,P15)
%
%  inputs:      fg       - function evaluation call
%               x0       - current design variable vector
%               f0       - objective value at x0
%               g0       - constraint values at x0
%               xmin     - minimum finite difference step
%               xmax     - maximum finite difference step
%               sx       - inline function to re-shape x
%               active   - active constraint index/boolean vector
%               Pn       - optional variables directly passed to fcnstr
%
%  outputs:     df       - finite difference objective gradient vector
%               dg       - finite difference constraint gradients matrix
%
%  Written by:   Robert A. Canfield
%
%  Created:      4/14/06
%  Modified:      5/5/08
%                9/26/15 - active constraints & complex step
%                10/2/16 - handle empty FinDiffRelStep

% Local variables
%
% dx....... Finite difference step
% f........ Perturbed objective value
% g........ Perturbed constraint values
% i........ Loop variable for current design variable perturbation

%--BEGIN
%
delx=1e-8;
if nargin<5 
   xmin=1e-8;
   xmax=1e-1;
else
   xmin=options.DiffMinChange;
   xmax=options.DiffMaxChange;
   if isfield(options,'FinDiffRelStep') && ~isempty(options.FinDiffRelStep)
      delx=options.FinDiffRelStep;
   end
end
if nargin<6 || isempty(active), active=1:numel(g0); end
% Less stringent relative change in x than 1.e-8 may be 
% needed for implicit functions that require numeric evaluation.
nac = numel(g0(active)); % number of constraints
dx = min( max(delx*abs(x0(:)),xmin), xmax );
dg = zeros(length(dx),nac);
df = zeros(size(dg,1),1);

% Forward Finite Difference or Complex Step loop.
for n=1:length(x0(:))
   x = x0;
   if strcmpi(options.FinDiffType,'complex')
      x(n) = x(n) + 1i*eps;
      [f,g] = fg(x,active);
      df(n) = imag(f) / eps;
      if nac
         dg(n,:) = imag(g(:)) / eps;
      end
   else
      x(n) = x(n) + dx(n);
      [f,g] = fg(x,active);
      df(n) = (f - f0) / dx(n);
      if nac
         dg(n,:) = (g(active) - g0(active)).'/ dx(n);
      end;
   end
end
end



%% Sub-function to transform user's fmincon function evauluation functions
function [f,g,nec] = fun2(x,Obj,Con,A,b,Aeq,beq)
f   = Obj(x);
if isempty(Con)
   g = [];
   h = [];
else
   [g,h] = Con(x);
end
if nargin<4
   A = [];
   b = [];
end
if nargin<6
   Aeq = [];
   beq = [];
end
if ~isempty(A) && ~isempty(b)
   c = A*x(:) - b(:);
   g = [g,c];
end
if ~isempty(Aeq) && ~isempty(beq)
   ceq = Aeq*x(:) - beq(:);
   h   = [h,ceq];
end
nec = length(h);
g = [h; g];
end



%% Sub-function to transform user's fmincon gradient evauluation functions
function [gradf,gradg] = grd2(x,Obj,Con,A,Aeq)
if isempty(Con)
   gradg = [];
   gradh = [];
else
   [~,gradf] = Obj(x);
   [~,~,gradg,gradh] = Con(x);
end
if nargin<4
   A = [];
end
if nargin<5
   Aeq = [];
end
if ~isempty(A)
   gradc = A;
   gradg = [gradg,gradc];
end
if ~isempty(Aeq)
   gradceq = Aeq;
   gradh   = [gradh,gradceq];
end
gradg = [gradh, gradg];
end