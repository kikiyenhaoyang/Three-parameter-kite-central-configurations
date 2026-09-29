function summary = verify_kite_regularity()
%VERIFY_KITE_REGULARITY Exact local regularity and uniform contraction.
%
% Requires Symbolic Math Toolbox. Run: verify_kite_regularity
% This certifies a closed box, not just a Jacobian at a single point.
% All interval endpoints and derivative bounds are exact symbolic rationals.
% This is a local mass-order-boundary result, not dynamical nondegeneracy.
if isempty(ver('symbolic'))
    error('verify_kite_regularity:MissingSymbolicToolbox', ...
        'This verification script requires Symbolic Math Toolbox.');
end
box={iv('0.50316','0.50317'),iv('0.38998','0.39000'),iv('0.199999','0.200001')};
x=box{1}; y=box{2}; z=box{3}; center=cell(1,3);
for k=1:3, center{k}=(box{k}.lo+box{k}.hi)/2; end
require([0<x.lo;x.lo<x.hi;x.hi<1;0<y.lo;y.lo<y.hi;y.hi<1; ...
    0<z.lo;z.lo<z.hi;z.hi<x.lo],'Basic geometry');
px=add(subtract(multiply(x,x),multiply(4,x)),1);
py=add(subtract(multiply(y,y),multiply(4,y)),1);
branch=multiply(multiply(3,x),x);
require([px.hi<0;py.hi<0],'Algebraic lower bounds');
require(branch.hi<1,'Lower admissibility branch');
F=subtract(multiply(subtract(1,multiply(y,y)),add(1,multiply(x,z))), ...
    multiply(multiply(2,y),subtract(x,z)));
G=subtract(multiply(x,subtract(1,multiply(z,z))),multiply(2,z));
Jgeom=subtract(add(multiply(subtract(1,multiply(x,x)),subtract(1,multiply(z,z))), ...
    multiply(multiply(4,x),z)), ...
    multiply(multiply(multiply(2,y),subtract(x,z)),add(1,multiply(x,z))));
require([F.lo>0;G.lo>0;Jgeom.lo>0],'Strict admissibility');
fprintf('PASS: the entire closed box B is strictly inside Omega_adm^J.\n');

[determinant,~]=kite_data(box{:});
strict_enclosure(determinant,-17,-16,'det(P)');
values=kite_masses(box); labels=[1 3 5];
fprintf('det(P)(B) subset %s\n',show(determinant));
for k=1:3
    strict_enclosure(values{k},'0.45','0.46',['m' num2str(labels(k))]);
    fprintf('m%d(B) subset %s\n',labels(k),show(values{k}));
end
fprintf('PASS: -17 < det(P) < -16 and 0.45 < m1,m3,m5 < 0.46 on B.\n');
DM=kite_jacobian(box); detDM=det3_interval(DM);
strict_enclosure(detDM,249,351,'det(DM)');
Df=cell(2,3);
for j=1:3
    Df{1,j}=subtract(DM{1,j},DM{2,j});
    Df{2,j}=subtract(DM{2,j},DM{3,j});
end
minor=subtract(multiply(Df{1,1},Df{2,2}),multiply(Df{1,2},Df{2,1}));
strict_enclosure(minor,115,122,'xy minor');
require(Df{1,1}.lo>11,'d_x(m1-m3)>11');
require(Df{2,2}.lo>8,'d_y(m3-m5)>8');
dx15=subtract(DM{1,1},DM{3,1});
require(dx15.lo>17,'d_x(m1-m5)>17');
fprintf('det(DM)(B) subset %s\n',show(detDM));
fprintf('det D_(x,y)(m1-m3,m3-m5)(B) subset %s\n',show(minor));
fprintf('PASS: 249 < det(DM) < 351 and 115 < xy minor < 122 on B.\n');
fprintf('PASS: d_x(m1-m3)>11, d_y(m3-m5)>8, d_x(m1-m5)>17 on B.\n');

% Rational preconditioner, uniform for every z in B_z.
preconditioner={iv(sym(73200)/10^6),iv(sym(23227)/10^6); ...
    iv(sym(-45490)/10^6),iv(sym(101217)/10^6)};
detR=subtract(multiply(preconditioner{1,1},preconditioner{2,2}), ...
    multiply(preconditioner{1,2},preconditioner{2,1}));
require([detR.lo==detR.hi;detR.lo~=0],'R must be invertible.');
E=cell(2,2); rowBounds=sym(zeros(1,2));
for i=1:2
    for j=1:2
        product=iv(0);
        for k=1:2, product=add(product,multiply(preconditioner{i,k},Df{k,j})); end
        E{i,j}=subtract(iv(double(i==j)),product);
    end
    for j=1:2, rowBounds(i)=rowBounds(i)+max([abs(E{i,j}.lo),abs(E{i,j}.hi)]); end
end
require(rowBounds<sym(3)/100,'Contraction bound');
mc=kite_masses(cellfun(@iv,center,'UniformOutput',false));
fc={subtract(mc{1},mc{2}),subtract(mc{2},mc{3})};
DMz=kite_jacobian({iv(center{1}),iv(center{2}),box{3}});
fz={subtract(DMz{1,3},DMz{2,3}),subtract(DMz{2,3},DMz{3,3})};
kCenter=matrix_vector(preconditioner,fc); kZ=matrix_vector(preconditioner,fz);
kError=matrix_vector(E,{subtract(box{1},center{1}),subtract(box{2},center{2})});
K=cell(1,2);
for i=1:2
    K{i}=add(subtract(subtract(iv(center{i}),kCenter{i}), ...
        multiply(kZ{i},subtract(box{3},center{3}))),kError{i});
    require([box{i}.lo<K{i}.lo;K{i}.lo<=K{i}.hi;K{i}.hi<box{i}.hi], ...
        'Uniform contraction image must lie strictly inside X.');
end
strict_enclosure(K{1},'0.5031634','0.5031677','K_x');
strict_enclosure(K{2},'0.3899882','0.3899906','K_y');
fprintf('K_x subset %s\nK_y subset %s\n',show(K{1}),show(K{2}));
fprintf('PASS: ||I-R D_(x,y)f||_infinity < 3/100 throughout B.\n');
fprintf('PASS: K is strictly inside B_x x B_y, uniformly for all z in B_z.\n');
fprintf('THEOREM: for every z in [0.199999,0.200001], exactly one (x,y)\n');
fprintf('in B_x x B_y has m1=m3=m5; these points form an analytic arc.\n');
fprintf('THEOREM: each equality surface is regular in int(B), and the\n');
fprintf('three surfaces meet pairwise transversely along this arc.\n');
fprintf('COROLLARY: every point of the arc in int(B) is approached by\n');
fprintf('all six strict positive axial mass orders.\n');
fprintf('SCOPE: this is a local kite result, not a global classification.\n');
summary=struct('detP',show(determinant),'detDM',show(detDM),'xy_minor',show(minor), ...
    'K_x',show(K{1}),'K_y',show(K{2}),'masses', ...
    {cellfun(@(v)show(v),values,'UniformOutput',false)});
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

function values = kite_masses(point)
[determinant,numerators]=kite_data(point{:});
values=cell(1,3);
for k=1:3, values{k}=divide(numerators{k},determinant); end
end

function J = kite_jacobian(point)
J=cell(3,3);
for axis=1:3
    variables=cell(1,3);
    for k=1:3, variables{k}=dv(point{k},double(k==axis)); end
    values=kite_masses(variables);
    for k=1:3, J{k,axis}=values{k}.derivative; end
end
end

function result = matrix_vector(a,b)
require(size(a,2)==numel(b),'Incompatible exact matrix/vector sizes.');
result=cell(1,size(a,1));
for i=1:size(a,1)
    result{i}=iv(0);
    for k=1:size(a,2), result{i}=add(result{i},multiply(a{i,k},b{k})); end
end
end

function [determinant,numerators] = kite_data(x,y,z)
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
determinant=multiply(multiply(alpha,subtract(1,alpha)),subtract(multiply(multiply(d1,d4),d5), ...
    multiply(multiply(d2,d3),d6)));
numerators={add(subtract(multiply(multiply(multiply(subtract(1,alpha),d1),d4),W), ...
    multiply(multiply(multiply(multiply(alpha,subtract(1,alpha)),d2),d6),V)), ...
    multiply(multiply(multiply(powerExact(subtract(1,alpha),2),d4),d6),U)), ...
    subtract(add(multiply(multiply(multiply(multiply(alpha,subtract(1,alpha)),d4),d5),U), ...
    multiply(multiply(multiply(alpha,d2),d3),W)),multiply(multiply(multiply(powerExact(alpha,2),d2), ...
    d5),V)),subtract(add(multiply(multiply(multiply(alpha,d1),d5),V), ...
    multiply(multiply(multiply(subtract(alpha,1),d3),d6),U)),multiply(multiply(d1,d3),W))};
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
