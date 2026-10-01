function summary = verify_kite_noninjectivity()
%VERIFY_KITE_NONINJECTIVITY Two distinct nonsingular kites for the same masses.
%
% Run: verify_kite_noninjectivity
% Or:  run_all('Only',{'N1'})
% Requires Symbolic Math Toolbox. All helpers are local to this file.
%
% The exact five-body mass vector is (1/50,1,1/100,1,3/10).
% A rational interval contraction proves a unique root of P*mu-b in each
% of two disjoint boxes. Both boxes are in Omega_CC, and both roots have
% det(DM) nonzero. Disjoint alpha=r15/r13 intervals rule out similarities
% and permutations preserving the masses. Thus M is not globally injective.
% No numerical root finder, floating-point tolerance or external data is used.
% The statement concerns P and DM, not dynamical/Hessian nondegeneracy.
if isempty(ver('symbolic'))
    error('verify_kite_noninjectivity:MissingSymbolicToolbox', ...
        'This verification script requires Symbolic Math Toolbox.');
end
verify_cartesian_identities();
target={'1/50','1/100','3/10'};
mu=cellfun(@exact,target,'UniformOutput',false);
require([mu{:}]>0,'The three axial masses must be strictly positive.');
centers={'0.424030037799','0.293530589884','0.141339619965'; ...
         '0.978349754805','0.295667726655','0.556483870165'};
radius=exact('0.000000001');
rawR=cell(1,2);
rawR{1}={ ...
    '0.001047474671','-0.026424332166','-0.056058901980'; ...
    '0.037648177754','-0.000370632171','-0.007221054729'; ...
    '0.000381464973','-0.021922999241','-0.013944812758'};
rawR{2}={ ...
    '0.002029398108','0.254893890548','0.034681863910'; ...
    '0.096746071781','-0.000698387694','-0.018196348006'; ...
    '-0.000351363946','0.169431240815','0.140372440606'};
summary=struct();
summary.mass_vector={'1/50','1','1/100','1','3/10'};
summary.axial_mass_vector=target;
summary.radius='1/1000000000';
alphaBoxes=cell(1,2); boxes=cell(1,2);
detPBounds={'-109','-108';'-0.00022','-0.00021'};
detDfBounds={'30959','30960';'-347','-346'};
detDMBounds={'284','285';'-1630000','-1620000'};
alphaBounds={'1.3067','1.3068';'2.0682','2.0684'};
for index=1:2
    center=cellfun(@exact,centers(index,:),'UniformOutput',false);
    box=cell(1,3);
    for j=1:3, box{j}=iv(center{j}-radius,center{j}+radius); end
    boxes{index}=box;
    [F,Jgeom,branchMargin]=verify_domain(box,index);
    [P,b,g]=kite_system(box{:});
    require([g.A.lo>0;g.B.lo>0;g.alpha.lo>1],'Positive shape denominators');
    detP=det3_interval(P);
    strict_enclosure(detP,detPBounds{index,1},detPBounds{index,2},'det(P)');
    masses=cell(1,3);
    for j=1:3
        replaced=P; replaced(:,j)=b(:);
        masses{j}=divide(det3_interval(replaced),detP);
        require(masses{j}.lo>0,'The entire box must have positive masses.');
    end

    Df=residual_jacobian(box,mu);
    detDf=det3_interval(Df);
    strict_enclosure(detDf,detDfBounds{index,1},detDfBounds{index,2},'det(Df)');
    R=cellfun(@exact,rawR{index},'UniformOutput',false);
    require(det(reshape([R{:}],3,3))~=0,'The rational preconditioner is invertible.');
    E=cell(3,3); rowBounds=sym(zeros(1,3));
    for i=1:3
        for j=1:3
            product=iv(0);
            for k=1:3, product=add(product,multiply(R{i,k},Df{k,j})); end
            E{i,j}=subtract(double(i==j),product);
            rowBounds(i)=rowBounds(i)+max([abs(E{i,j}.lo),abs(E{i,j}.hi)]);
        end
    end
    q=max(rowBounds);
    require(q<sym(1)/100000,'The contraction factor must be below 1e-5.');
    fc=fixed_mass_residual(center,mu);
    correction=matrix_vector(R,fc);
    errorBox=matrix_vector(E,{iv(-radius,radius),iv(-radius,radius),iv(-radius,radius)});
    image=cell(1,3);
    for i=1:3
        image{i}=add(subtract(center{i},correction{i}),errorBox{i});
        require([box{i}.lo<image{i}.lo;image{i}.lo<=image{i}.hi;image{i}.hi<box{i}.hi], ...
            'The contraction image must lie strictly inside the box.');
        require([center{i}-radius/100<image{i}.lo;image{i}.hi<center{i}+radius/100], ...
            'The image must lie within 1e-11 of the rational center.');
    end
    % Banach's theorem now gives exactly one p in the interior with P(p)mu=b(p).
    % At that root only: Df=-P*DM, hence det(DM)=-det(Df)/det(P).
    detDMAtRoot=divide(negate(detDf),detP);
    strict_enclosure(detDMAtRoot,detDMBounds{index,1},detDMBounds{index,2}, ...
        'det(DM) at the certified root');
    alphaBoxes{index}=g.alpha;
    strict_enclosure(g.alpha,alphaBounds{index,1},alphaBounds{index,2},'alpha');
    lambda=add(add(add(multiply(mu{1},powerExact(g.r12,-3)), ...
        multiply(mu{2},powerExact(g.r23,-3))),multiply(mu{3},powerExact(g.r25,-3))), ...
        multiply(2,powerExact(g.r24,-3)));
    require(lambda.lo>0,'The central-configuration multiplier is positive.');
    fprintf('PASS box %d: a unique exact root with masses (1/50,1,1/100,1,3/10).\n',index);
    fprintf('det(P): %s; det(Df): %s\n',show(detP,15),show(detDf,15));
    fprintf('det(DM) at the root: %s\n',show(detDMAtRoot,15));
    fprintf('Contraction factor <= %s < 1e-5\n',decimal(q,true,15));
    fprintf('alpha: %s\n',show(g.alpha,15));
    for j=1:3, fprintf('Root coordinate %d in %s\n',j,show(image{j},15)); end
    summary.boxes(index)=struct('center',{centers(index,:)}, ...
        'box',{cellfun(@(v)show(v,12),box,'UniformOutput',false)}, ...
        'preconditioner',{rawR{index}},'detP',show(detP,15),'detDf',show(detDf,15), ...
        'detDM_at_root',show(detDMAtRoot,15),'contraction_bound',decimal(q,true,15), ...
        'image',{cellfun(@(v)show(v,15),image,'UniformOutput',false)}, ...
        'alpha',show(g.alpha,15), ...
        'masses_on_box',{cellfun(@(v)show(v,15),masses,'UniformOutput',false)}, ...
        'F',show(F,15),'J',show(Jgeom,15),'branch_margin',show(branchMargin,15), ...
        'lambda_enclosure',show(lambda,15));
end
firstBox=boxes{1};secondBox=boxes{2};
require(firstBox{1}.hi<secondBox{1}.lo,'The two parameter boxes are disjoint.');
require(alphaBoxes{1}.hi<alphaBoxes{2}.lo,'The similarity invariant alpha separates the roots.');
% The masses of bodies 1,3,5 are pairwise distinct and differ from 1.
% Every mass-preserving permutation fixes them, so r15/r13 is invariant.
require([mu{1}~=mu{2};mu{1}~=mu{3};mu{2}~=mu{3}; ...
    mu{1}~=1;mu{2}~=1;mu{3}~=1],'Unique axial mass labels');
summary.alpha_intervals_disjoint=true;
summary.exact_roots_certified=2;
summary.global_noninjectivity=true;
summary.local_inverse_at_both_roots=true;
fprintf('THEOREM: the normalized nonsingular mass map is not globally injective.\n');
fprintf('The two mass-preserving similarity classes are distinct.\n');
fprintf('COROLLARY: every mass vector in some open neighborhood of the target has at least two shapes.\n');
end

function [F,Jgeom,margin] = verify_domain(box,index)
x=box{1};y=box{2};z=box{3};
require([0<x.lo;x.lo<x.hi;x.hi<1;0<y.lo;y.lo<y.hi;y.hi<1; ...
    0<z.lo;z.lo<z.hi;z.hi<x.lo],'Basic shape bounds');
px=add(subtract(multiply(x,x),multiply(4,x)),1);
py=add(subtract(multiply(y,y),multiply(4,y)),1);
require([px.hi<0;py.hi<0],'x,y exceed 2-sqrt(3)');
F=subtract(multiply(subtract(1,multiply(y,y)),add(1,multiply(x,z))), ...
    multiply(multiply(2,y),subtract(x,z)));
Jgeom=subtract(add(multiply(subtract(1,multiply(x,x)),subtract(1,multiply(z,z))), ...
    multiply(multiply(4,x),z)), ...
    multiply(multiply(multiply(2,y),subtract(x,z)),add(1,multiply(x,z))));
require([F.lo>0;Jgeom.lo>0],'F and J are strictly positive');
branch=multiply(multiply(3,x),x);
if index==1
    require(branch.hi<1,'First box is in the lower branch');
    margin=subtract(multiply(x,subtract(1,multiply(z,z))),multiply(2,z));
else
    require(branch.lo>1,'Second box is in the upper branch');
    h=divide(subtract(x,z),add(1,multiply(x,z)));
    require([0<h.lo;h.hi<1],'Half-angle difference bounds');
    margin=negate(add(subtract(multiply(h,h),multiply(4,h)),1));
end
require(margin.lo>0,'Strict branch admissibility');
end

function values = fixed_mass_residual(point,mu)
[P,b]=kite_system(point{:});
values=cell(1,3);
for i=1:3
    values{i}=negate(b{i});
    for j=1:3, values{i}=add(values{i},multiply(P{i,j},mu{j})); end
end
end

function J = residual_jacobian(point,mu)
J=cell(3,3);
for axis=1:3
    variables=cell(1,3);
    for j=1:3, variables{j}=dv(point{j},double(j==axis)); end
    values=fixed_mass_residual(variables,mu);
    for i=1:3, J{i,axis}=values{i}.derivative; end
end
end

function verify_cartesian_identities()
% Exact polynomial identities after clearing denominators. A and B are
% certified strictly positive on each interval box in the main verifier.
% Free symbols must not pass through the guarded interval division helpers.
x=sym('x','real');
y=sym('y','real');
z=sym('z','real');
A=(x+y)*(1-x*y);
B=(x-z)*(1+x*z);
R=y*z*(1+x*x)^2;
C=A*B+R;
aNumerator=x*(1-y*y);
hNumerator=2*x*y;
r12Numerator=x*(1+y*y);
r23Numerator=y*(1+x*x);
r25Numerator=x*y*(1+x*x)*(1+z*z);
equalZero([2*y*(1-x*x)+2*aNumerator-2*A; ...
    C-A*B-R; ...
    r12Numerator^2-aNumerator^2-hNumerator^2; ...
    r23Numerator^2-(A-aNumerator)^2-hNumerator^2; ...
    r25Numerator^2-(C-aNumerator*B)^2-(hNumerator*B)^2], ...
    'Denominator-cleared Cartesian distance identities');
a=sym('a','real');
alpha=sym('alpha','real');
m1=sym('m1','real');
m3=sym('m3','real');
m5=sym('m5','real');
s12=sym('s12','real');
s23=sym('s23','real');
s24=sym('s24','real');
s25=sym('s25','real');
s35=sym('s35','real');
s15=sym('s15','real');
lambda=m1*s12+m3*s23+m5*s25+2*s24;
a2=-m1*a*s12+m3*(1-a)*s23+m5*(alpha-a)*s25;
a1=2*a*s12+m3+m5*alpha*s15;
a3=-m1+2*(a-1)*s23+m5*(alpha-1)*s35;
a5=-m1*alpha*s15-m3*(alpha-1)*s35+2*(a-alpha)*s25;
f1=(1-s23)*m3+alpha*(s15-s25)*m5-2*a*(s24-s12);
f2=(s12-1)*m1+(1-alpha)*(s25-s35)*m5-2*(1-a)*(s23-s24);
f3=alpha*(s12-s15)*m1+(alpha-1)*(s23-s35)*m3-2*(alpha-a)*(s25-s24);
equalZero([a1-a2-lambda*a-f1;a3-a2+lambda*(1-a)-f2; ...
    a5-a2+lambda*(alpha-a)-f3],'Full Cartesian reduction identities');
% For positive masses and distances, lambda>0. Reflection gives the y
% equations. When f=0, all a_i+lambda*q_i coincide; total force balance
% then identifies their common value as lambda times the center of mass.
fprintf('PASS: symbolic distance and full Cartesian reduction identities.\n');
end

function [P,b,g] = kite_system(x,y,z)
% Exact interval/dual arithmetic in the original evaluation order.
A=multiply(add(x,y),subtract(1,multiply(x,y)));
B=multiply(subtract(x,z),add(1,multiply(x,z)));
R=multiply(multiply(y,z),powerExact(add(1,multiply(x,x)),2));
alpha=add(1,divide(R,multiply(A,B)));
beta=divide(multiply(multiply(2,y),subtract(1,multiply(x,x))),A);
gamma=divide(multiply(multiply(2,x),subtract(1,multiply(y,y))),A);
r12=divide(multiply(x,add(1,multiply(y,y))),A);
r23=divide(multiply(y,add(1,multiply(x,x))),A);
r24=divide(multiply(multiply(4,x),y),A);
r25=divide(multiply(multiply(multiply(x,y),add(1,multiply(x,x))),add(1,multiply(z,z))),multiply(A,B));
r35=divide(R,multiply(A,B));
S12=powerExact(r12,-3);
S23=powerExact(r23,-3);
S24=powerExact(r24,-3);
S25=powerExact(r25,-3);
S35=powerExact(r35,-3);
S15=powerExact(alpha,-3);
d1=subtract(1,S23);
d2=subtract(S15,S25);
d3=subtract(S12,1);
d4=subtract(S25,S35);
d5=subtract(S12,S15);
d6=subtract(S23,S35);
U=multiply(subtract(S24,S12),gamma);
V=multiply(subtract(S23,S24),beta);
W=multiply(subtract(S25,S24),subtract(add(beta,multiply(2,alpha)),2));
P={0,d1,multiply(alpha,d2); ...
   d3,0,multiply(subtract(1,alpha),d4); ...
   multiply(alpha,d5),multiply(subtract(alpha,1),d6),0};
b={U;V;W};
g=struct('A',A,'B',B,'alpha',alpha,'beta',beta,'gamma',gamma, ...
    'r12',r12,'r23',r23,'r24',r24,'r25',r25,'r35',r35);
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

function value = iv(lo,hi)
% Closed interval with exact symbolic rational endpoints.
if isstruct(lo) && isfield(lo,'lo') && isfield(lo,'hi')
    value=lo;
    return
end
if nargin<2, hi=lo; end
lo=exact(lo); hi=exact(hi);
require(lo<=hi,'Reversed interval.');
value=struct('lo',lo,'hi',hi);
end

function c = add(a,b)
if isDual(a) || isDual(b)
    a=dv(a); b=dv(b);
    c=dv(add(a.value,b.value),add(a.derivative,b.derivative));
elseif isstruct(a) || isstruct(b)
    a=iv(a); b=iv(b);
    c=iv(a.lo+b.lo,a.hi+b.hi);
else
    c=exact(a)+exact(b);
end
end

function c = subtract(a,b)
c=add(a,negate(b));
end

function c = multiply(a,b)
if isDual(a) || isDual(b)
    a=dv(a); b=dv(b);
    c=dv(multiply(a.value,b.value), ...
        add(multiply(a.derivative,b.value),multiply(a.value,b.derivative)));
elseif isstruct(a) || isstruct(b)
    a=iv(a); b=iv(b);
    endpoints=[a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi];
    c=iv(min(endpoints),max(endpoints));
else
    c=exact(a)*exact(b);
end
end

function text = show(value,places)
if nargin<2, places=12; end
text=['[' decimal(value.lo,false,places) ', ' decimal(value.hi,true,places) ']'];
end

function strict_enclosure(value,lo,hi,label)
require([exact(lo)<value.lo; value.lo<=value.hi; value.hi<exact(hi)], ...
    [label ': enclosure failed.']);
end

function result = matrix_vector(a,b)
require(size(a,2)==numel(b),'Incompatible exact matrix/vector sizes.');
result=cell(1,size(a,1));
for i=1:size(a,1)
    result{i}=iv(0);
    for k=1:size(a,2), result{i}=add(result{i},multiply(a{i,k},b{k})); end
end
end

function result = det3_interval(a)
% Exact interval/dual arithmetic in the original evaluation order.
result=add(subtract(multiply(a{1,1},subtract(multiply(a{2,2},a{3,3}),multiply(a{2,3},a{3,2}))), ...
    multiply(a{1,2},subtract(multiply(a{2,1},a{3,3}),multiply(a{2,3},a{3,1})))),multiply(a{1,3}, ...
    subtract(multiply(a{2,1},a{3,2}),multiply(a{2,2},a{3,1}))));
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

function value = dv(a,derivative)
% Interval value and enclosure of one partial derivative.
if isDual(a), value=a; return; end
if nargin<2, derivative=0; end
value=struct('value',iv(a),'derivative',iv(derivative));
end

function yes = isDual(a)
yes=isstruct(a) && isfield(a,'value') && isfield(a,'derivative');
end

function c = negate(a)
if isDual(a)
    c=dv(negate(a.value),negate(a.derivative));
elseif isstruct(a)
    c=iv(-a.hi,-a.lo);
else
    c=-exact(a);
end
end

function text = decimal(q,upper,places)
% Directed decimal display formed by exact symbolic integer division.
if nargin<2, upper=false; end
if nargin<3, places=12; end
require(isnumeric(places) && isscalar(places) && places>=0 && fix(places)==places, ...
    'Invalid number of decimal places.');
q=exact(q); scale=sym(10)^places;
if upper, rounded=ceil(q*scale); else, rounded=floor(q*scale); end
negative=proved(rounded<0); digits=char(abs(rounded));
if places>0
    digits=[repmat('0',1,max(0,places+1-numel(digits))) digits];
    text=[digits(1:end-places) '.' digits(end-places+1:end)];
else
    text=digits;
end
if negative, text=['-' text]; end
end

function c = divide(a,b)
c=multiply(a,powerExact(b,-1));
end

function c = powerExact(a,p)
require(isnumeric(p) && isscalar(p) && isreal(p) && isfinite(p) && fix(p)==p, ...
    'Expected an integer power.');
if isDual(a)
    if p==0, c=dv(1); return; end
    c=dv(powerExact(a.value,p), ...
        multiply(multiply(p,powerExact(a.value,p-1)),a.derivative));
elseif isstruct(a)
    if p<0, c=powerExact(reciprocal(a),-p); return; end
    if p==0, c=iv(1); return; end
    if mod(p,2)==1, c=iv(a.lo^p,a.hi^p); return; end
    lower=min([a.lo^p,a.hi^p]);
    if proved((a.lo<=0) & (0<=a.hi)), lower=sym(0); end
    c=iv(lower,max([a.lo^p,a.hi^p]));
else
    a=exact(a);
    if p<0, require(a~=0,'Zero rational denominator.'); end
    c=a^p;
end
end

function yes = proved(condition)
if isa(condition,'sym')
    yes=all(isAlways(condition(:),'Unknown','false'));
else
    yes=all(logical(condition(:)));
end
end

function c = reciprocal(a)
% Called only on a plain interval or a scalar symbolic value.
if isstruct(a)
    require((a.hi<0) | (a.lo>0),'Interval denominator contains zero.');
    c=iv(1/a.hi,1/a.lo);
else
    a=exact(a); require(a~=0,'Zero rational denominator.');
    c=1/a;
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
