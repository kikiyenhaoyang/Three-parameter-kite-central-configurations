function summary = verify_region_intersections()
%VERIFY_REGION_INTERSECTIONS Exact positive-mass intersections in A,B,F,G.
% Each segment has exactly one H zero. No global component count is claimed.
% Requires Symbolic Math Toolbox. All helpers are local to this file.
% Comparisons use exact symbolic rationals; no floating-point tolerances.
if isempty(ver('symbolic'))
    error('verify_region_intersections:MissingSymbolicToolbox', ...
        'This verification script requires Symbolic Math Toolbox.');
end
segments={'A','0.31','0.4710','0.4711'; 'B','0.35','0.5878','0.5879'; ...
    'F','0.40','0.7275','0.7276'; 'G','0.45','0.7705','0.7706'};
orders={ {'15','25','12','24','35','13'}, {'15','25','12','35','24','13'}, ...
    {'15','25','35','12','24','13'}, {'15','25','35','12','13','24'} };
summary=struct();
for k=1:4
    region=segments{k,1}; point=rational_tuple(segments(k,2:4));
    u=point{1}; v=point{2}; w=point{3};
    require(0<u & u<v & v<w & w<1,'Basic domain inequalities failed.');
    require(u*u-4*u+1<0 & 3*u*u<1,'Algebraic angle bounds failed.');
    f=fields_uv(iv(u),iv(v,w));
    r=distances_uv(iv(u),iv(v,w)); order=orders{k};
    for j=1:5
        gap=subtract(r.(['r' order{j+1}]),r.(['r' order{j}]));
        require(gap.lo>exact(3,100),'The segment does not have the required strict distance order.');
    end
    require(0<f.t.lo & f.t.hi<exact(1,2),'Invalid axial position.');
    ca=subtract(f.C,f.A);
    require(f.A.lo>1 & ca.lo>0 & f.D.lo>0,'A required denominator sign failed.');
    require(f.m1.lo>exact(1,10) & f.m3.lo>exact(1,10) & f.m5.lo>exact(1,10), ...
        'The entire segment must have strictly positive masses.');
    left=fields_uv(u,v); right=fields_uv(u,w);
    require(left.H>exact(1,10000),'Lower endpoint H must exceed 1/10000.');
    require(right.H<-exact(1,10000),'Upper endpoint H must be below -1/10000.');
    differentiated=fields_uv(dv(iv(u)),dv(iv(v,w),1));
    hv=differentiated.H.derivative;
    require(hv.hi<-6,'The v derivative must be below -6 throughout the segment.');
    endpoints={left,right};
    for j=1:2
        e=endpoints{j};
        equalZero(e.m5-e.m1-e.H/e.D,'Mass-difference identity');
        t=e.t; b=e.b; A=e.A; B=e.B; C=e.C;
        m1=e.m1; m3=e.m3; m5=e.m5;
        residuals={ ...
            (1-A)*m3+t*(t^(-3)-C)*m5-(B-A), ...
            (A-1)*m1+b*(C-b^(-3))*m5-(A-B), ...
            t*(A-t^(-3))*m1+b*(b^(-3)-A)*m3-(B-C)*(1-2*t)};
        equalZero([residuals{:}],'Original reduced-system residuals');
    end
    fprintf('PASS %s: u=%s, %s <= v <= %s.\n',segments{k,:});
    fprintf('  H(lower) in %s\n',show(iv(left.H),10));
    fprintf('  H(upper) in %s\n',show(iv(right.H),10));
    fprintf('  H_v(segment) in %s\n',show(hv,10));
    fprintf('  m1(segment) in %s\n',show(f.m1,10));
    fprintf('  m3(segment) in %s\n',show(f.m3,10));
    fprintf('  m5(segment) in %s\n',show(f.m5,10));
    fprintf('  Every consecutive distance gap exceeds 3/100.\n');
    fprintf('  Thus this positive-mass segment contains exactly one zero of H.\n');
    summary.regions(k)=struct('region',region, ...
        'H_lower',show(iv(left.H),10),'H_upper',show(iv(right.H),10), ...
        'Hv',show(hv,10),'m1',show(f.m1,10),'m3',show(f.m3,10),'m5',show(f.m5,10));
end
fprintf('VERIFIED: H=0 meets the OPEN regions A, B, F and G in positive masses.\n');
fprintf('Each certified zero has H_v<-6 and hence lies on a local analytic arc.\n');
fprintf('The independent closure exclusion is checked separately by certificate R3.\n');
fprintf('No global component count or global graph property is asserted.\n');
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

function value = dv(a,derivative)
% Interval value and enclosure of one partial derivative.
if isDual(a), value=a; return; end
if nargin<2, derivative=0; end
value=struct('value',iv(a),'derivative',iv(derivative));
end

function c = subtract(a,b)
c=add(a,negate(b));
end

function text = show(value,places)
if nargin<2, places=12; end
text=['[' decimal(value.lo,false,places) ', ' decimal(value.hi,true,places) ']'];
end

function values = rational_tuple(raw)
if isstring(raw), raw=cellstr(raw); end
require(iscell(raw),'Rational coordinates must be JSON strings.');
values=cell(1,numel(raw));
for k=1:numel(raw)
    require(ischar(raw{k}) || (isstring(raw{k}) && isscalar(raw{k})), ...
        'Rational coordinates must be strings.');
    values{k}=exact(raw{k});
end
end

function result = fields_uv(u,v)
% Exact interval/dual arithmetic in the original evaluation order.
t=divide(multiply(subtract(v,u),add(1,multiply(u,v))),multiply(multiply(2,v),subtract(1,multiply(u,u))));
b=subtract(1,t);
q=multiply(t,b);
A=divide(multiply(8,powerExact(subtract(1,multiply(u,u)),3)),powerExact(add(1,multiply(u,u)),3));
B=divide(powerExact(subtract(1,multiply(u,u)),3),multiply(8,powerExact(u,3)));
C=divide(multiply(multiply(8,powerExact(v,3)),powerExact(subtract(1,multiply(u,u)),3)), ...
    multiply(powerExact(u,3),powerExact(add(1,multiply(v,v)),3)));
N=add(subtract(A,B),multiply(multiply(q,q),subtract(subtract(multiply(A,A),B),multiply(subtract(A,1),C))));
D=multiply(multiply(multiply(multiply(b,b),subtract(A,1)),subtract(C,A)),subtract(1,q));
F1=subtract(multiply(multiply(multiply(multiply(b,b),subtract(A,B)),subtract(C,A)),subtract(1,q)), ...
    multiply(subtract(multiply(powerExact(b,3),C),1),N));
F5=multiply(multiply(multiply(b,b),subtract(A,1)),N);
H=subtract(F5,F1);
m1=divide(F1,D);
m5=divide(F5,D);
m3=add(m1,multiply(divide(multiply(subtract(1,multiply(2,t)),add(powerExact(q,-2),C)),subtract(A,1)),m5));
result=struct('t',t,'b',b,'q',q,'A',A,'B',B,'C',C,'N',N,'D',D,'F1',F1,'F5',F5,'H',H,'m1',m1,'m3',m3, ...
    'm5',m5);
end

function result = distances_uv(u,v)
% Exact interval/dual arithmetic in the original evaluation order.
t=divide(multiply(subtract(v,u),add(1,multiply(u,v))),multiply(multiply(2,v),subtract(1,multiply(u,u))));
result=struct('r15',t,'r35',subtract(1,t),'r13',add(multiply(u,0),1),'r12',divide(add(1,multiply(u, ...
    u)),multiply(2,subtract(1,multiply(u,u)))),'r24',divide(multiply(2,u),subtract(1,multiply(u, ...
    u))),'r25',divide(multiply(u,add(1,multiply(v,v))),multiply(multiply(2,v),subtract(1,multiply(u,u)))));
end

function yes = isDual(a)
yes=isstruct(a) && isfield(a,'value') && isfield(a,'derivative');
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

function c = divide(a,b)
c=multiply(a,powerExact(b,-1));
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
