function summary = verify_kite_mass_orders(varargin)
%VERIFY_KITE_MASS_ORDERS Exact certificates for all six axial mass orders.
%
% Requires Symbolic Math Toolbox. Run from MATLAB with:
%     verify_kite_mass_orders
%     verify_kite_mass_orders('Exact',true)
%
% Main function plus local helper functions, following the supplied example.
% All inequalities and residuals use exact symbolic rational arithmetic.
if isempty(ver('symbolic'))
    error('verify_kite_mass_orders:MissingSymbolicToolbox', ...
        'This verification script requires Symbolic Math Toolbox.');
end
p=inputParser;
addParameter(p,'Exact',false,@(v)islogical(v)&&isscalar(v));
parse(p,varargin{:});

syms x y z real
vars=[x,y,z];
[P,b,geometry]=kite_system(x,y,z);
orders=[1 3 5;1 5 3;3 1 5;3 5 1;5 1 3;5 3 1];
points={'11/20','2/5','1/4'; '7/20','9/20','3/20'; ...
    '7/10','2/5','3/10'; '7/10','1/2','1/10'; ...
    '1/2','1/2','1/5'; '3/5','2/5','1/5'};
bounds={ [22 23;348 349;927 928], [69 70;1609 1610;454 455], ...
    [494 495;315 316;642 643], [103076 103077;964 965;1522 1523], ...
    [987 988;1669 1670;263 264], [2119 2120;520 521;177 178] };
require(isequal(sortrows(orders),sortrows(perms([1 3 5]))), ...
    'The witnesses do not cover all six strict mass orders.');
summary=struct();
for k=1:6
    point=cellfun(@exact,points(k,:),'UniformOutput',false);
    p0=[point{:}];
    x0=p0(1); y0=p0(2); z0=p0(3);
    require([0<x0; x0<1; 0<y0; y0<1; 0<z0; z0<x0], ...
        'Basic shape inequalities failed.');
    require([x0*x0-4*x0+1<0; y0*y0-4*y0+1<0], ...
        'Algebraic angle lower bounds failed.');
    g=struct(); names=fieldnames(geometry);
    for j=1:numel(names), g.(names{j})=subs(geometry.(names{j}),vars,p0); end
    require([g.F>0;g.J>0],'F or J is not strictly positive.');
    require(3*x0*x0~=1,'Witness is on the branch interface.');
    if proved(3*x0*x0<1)
        require(g.G>0,'G is not strictly positive.');
    else
        h=(x0-z0)/(1+x0*z0);
        require([0<h;h<1;h*h-4*h+1<0],'H is not strictly positive.');
    end
    require([g.A>0;g.B>0;g.alpha>1],'Invalid denominators or alpha.');
    equalZero(g.beta+g.gamma-2,'Axial-coordinate identity');

    P0=subs(P,vars,p0); b0=subs(b,vars,p0);
    % Keep an exact rational scalar; factor(rational) can return a vector.
    detP=simplifyFraction(det(P0));
    require(isscalar(detP),'det P must be scalar.');
    require(detP<0,'The certified determinant sign failed.');
    equalZero(detP-g.alpha*(1-g.alpha)* ...
        (g.d1*g.d4*g.d5-g.d2*g.d3*g.d6),'Determinant formula');
    m0=simplify(P0\b0);
    equalZero(P0*m0-b0,'Exact linear-system residual');

    % Independently retain both original Cramer checks.
    cramer=sym(zeros(3,1));
    for j=1:3
        replaced=P0; replaced(:,j)=b0;
        cramer(j)=simplifyFraction(det(replaced)/detP);
    end
    % Assign scalar entries explicitly: whitespace before unary +/- inside
    % a bracket literal can otherwise create extra columns in MATLAB.
    compact=sym(zeros(3,1));
    compact(1)=(1-g.alpha)*g.d1*g.d4*g.W-g.alpha*(1-g.alpha)*g.d2*g.d6*g.V + ...
        (1-g.alpha)^2*g.d4*g.d6*g.U;
    compact(2)=g.alpha*(1-g.alpha)*g.d4*g.d5*g.U+g.alpha*g.d2*g.d3*g.W - ...
        g.alpha^2*g.d2*g.d5*g.V;
    compact(3)=g.alpha*g.d1*g.d5*g.V+(g.alpha-1)*g.d3*g.d6*g.U-g.d1*g.d3*g.W;
    equalZero(m0-cramer,'Cramer solution');
    equalZero(m0-compact/detP,'Compact Cramer formulas');
    limits=bounds{k};
    for j=1:3
        require([sym(limits(j,1))<1000*m0(j);1000*m0(j)<sym(limits(j,2))], ...
            'Exact positive-mass table bounds');
    end
    ix=(orders(k,:)+1)/2;
    require([0<m0(ix(1));m0(ix(1))<m0(ix(2));m0(ix(2))<m0(ix(3))], ...
        'The prescribed strict positive mass order failed.');
    multiplier=verify_cartesian(g,m0);

    fprintf('PASS m%d < m%d < m%d; (x,y,z)=(%s); 1000*(m1,m3,m5) in %s\n', ...
        orders(k,1),orders(k,2),orders(k,3),strjoin(points(k,:),','),mat2str(limits));
    strings=cell(1,3);
    for j=1:3, strings{j}=rationalText(m0(j)); end
    if p.Results.Exact
        fprintf('Mass vector, exact =\n'); disp(m0.');
        fprintf('Mass vector, decimals (display only) =\n'); disp(vpa(m0.',15));
        fprintf('det P = %s\nlambda = %s\n',rationalText(detP),rationalText(multiplier));
    end
    summary.witnesses(k)=struct('order',orders(k,:),'point',{points(k,:)}, ...
        'masses',{strings},'determinant',rationalText(detP),'multiplier',rationalText(multiplier));
end
summary.witness_count=6; summary.cartesian_residual_count=60;
fprintf('PASS: all six rational witnesses are strictly admissible and nonsingular.\n');
fprintf('PASS: all mass bounds and all 60 Cartesian CC residuals are exact.\n');
fprintf('The nonempty open neighborhoods follow by continuity, as in the proposition.\n');
end

function [P,b,g] = kite_system(x,y,z)
A=(x+y)*(1-x*y);
B=(x-z)*(1+x*z);
R=y*z*(1+x*x)^2;
alpha=1+R/(A*B);
beta=2*y*(1-x*x)/A;
gamma=2*x*(1-y*y)/A;
r12=x*(1+y*y)/A;
r23=y*(1+x*x)/A;
r24=4*x*y/A;
r25=x*y*(1+x*x)*(1+z*z)/(A*B);
r35=R/(A*B);
s12=r12^(-3); s23=r23^(-3); s24=r24^(-3);
s25=r25^(-3); s35=r35^(-3); s15=alpha^(-3);
d1=1-s23; d2=s15-s25; d3=s12-1;
d4=s25-s35; d5=s12-s15; d6=s23-s35;
U=(s24-s12)*gamma;
V=(s23-s24)*beta;
W=(s25-s24)*(beta+2*alpha-2);
P=[0,d1,alpha*d2;d3,0,(1-alpha)*d4;alpha*d5,(alpha-1)*d6,0];
b=[U;V;W];
F=(1-y*y)*(1+x*z)-2*y*(x-z);
G=x*(1-z*z)-2*z;
J=(1-x*x)*(1-z*z)+4*x*z-2*y*(x-z)*(1+x*z);
g=struct('A',A,'B',B,'F',F,'G',G,'J',J,'alpha',alpha,'beta',beta,'gamma',gamma, ...
    'r12',r12,'r23',r23,'r24',r24,'r25',r25,'r35',r35, ...
    'd1',d1,'d2',d2,'d3',d3,'d4',d4,'d5',d5,'d6',d6,'U',U,'V',V,'W',W);
end

function multiplier = verify_cartesian(g,m0)
positions=[0,0;g.gamma/2,g.r24/2;1,0;g.gamma/2,-g.r24/2;g.alpha,0];
masses=[m0(1);sym(1);m0(2);sym(1);m0(3)];
pairs=[1 2;1 3;1 4;1 5;2 3;2 4;2 5;3 4;3 5;4 5];
lengths=[g.r12,sym(1),g.r12,g.alpha,g.r23,g.r24,g.r25,g.r23,g.r35,g.r25];
distances=sym(zeros(5));
for e=1:10
    i=pairs(e,1); j=pairs(e,2); distance=lengths(e);
    require(distance>0,'A collision or negative distance occurred.');
    equalZero(sum((positions(i,:)-positions(j,:)).^2)-distance^2, ...
        'Exact Cartesian distance');
    distances(i,j)=distance; distances(j,i)=distance;
end
center=(masses.'*positions)/sum(masses);
acceleration=sym(zeros(5,2));
for i=1:5
    for j=1:5
        if j~=i
            acceleration(i,:)=acceleration(i,:)+ ...
                masses(j)*(positions(j,:)-positions(i,:))/distances(i,j)^3;
        end
    end
end
multiplier=simplifyFraction(acceleration(1,1)-acceleration(3,1));
require(multiplier>0,'The central-configuration multiplier is not positive.');
for i=1:5
    equalZero(acceleration(i,:)+multiplier*(positions(i,:)-center), ...
        'Original Cartesian central-configuration residual');
end
end

function require(condition,label)
% Accept only statements proved true; undecidable statements fail closed.
if isa(condition,'sym')
    truth=false(size(condition));
    for k=1:numel(condition)
        truth(k)=isAlways(condition(k),'Unknown','false');
    end
else
    truth=logical(condition);
end
if ~all(truth(:))
    error('kite:CheckFailed','%s',label);
end
end

function yes = proved(condition)
if isa(condition,'sym')
    yes=all(isAlways(condition(:),'Unknown','false'));
else
    yes=all(logical(condition(:)));
end
end

function equalZero(expression,label)
% Prove each rational expression is zero via its numerator.
for k=1:numel(expression)
    reduced=simplifyFraction(expression(k));
    [numerator,~]=numden(reduced);
    numerator=simplify(expand(numerator));
    require(numerator==0,label);
end
end

function q = exact(a,b)
% Parse integers, fractions and finite decimal strings without double.
if nargin==2
    a=exact(a); b=exact(b);
    require(b~=0,'Zero rational denominator.');
    q=a/b;
    return
end
if isa(a,'sym')
    require(isscalar(a),'Expected a scalar exact value.');
    q=a;
    return
end
if isstring(a) && isscalar(a), a=char(a); end
if isnumeric(a) || islogical(a)
    require(isscalar(a) && isreal(a),'Expected a real scalar integer.');
    if isinteger(a)
        q=sym(char(string(a)));
    else
        require(isfinite(a) && fix(a)==a && abs(a)<=flintmax, ...
            'Use a string or sym(integer)/integer for fractional input.');
        q=sym(a);
    end
    return
end
require(ischar(a) && isrow(a),'Expected an exact rational string.');
a=strtrim(a); slash=strfind(a,'/');
if ~isempty(slash)
    require(numel(slash)==1,'Invalid rational string.');
    n=a(1:slash-1); d=a(slash+1:end);
    require(~isempty(regexp(n,'^[+-]?[0-9]+$','once')) && ...
        ~isempty(regexp(d,'^[+-]?[0-9]+$','once')),'Invalid fraction.');
    denominator=sym(d); require(denominator~=0,'Zero rational denominator.');
    q=sym(n)/denominator;
else
    require(~isempty(regexp(a,'^[+-]?[0-9]+(\.[0-9]+)?$','once')), ...
        'Expected an integer or finite decimal string.');
    dot=strfind(a,'.');
    if isempty(dot)
        q=sym(a);
    else
        digits=numel(a)-dot; a(dot)=[];
        q=sym(a)/sym(10)^digits;
    end
end
end

function text = rationalText(q)
[n,d]=numden(simplifyFraction(q));
text=char(n);
if ~proved(d==1), text=[text '/' char(d)]; end
end
