function summary = verify_exclusion(path)
%VERIFY_EXCLUSION Exact exclusion of relative closures C,D,E,H,I.
% Uses explicit H, geometric inequalities, and the exact corner factors.
% Requires Symbolic Math Toolbox. All helpers are local to this file.
% Comparisons use exact symbolic rationals; no floating-point tolerances.
if isempty(ver('symbolic'))
    error('verify_exclusion:MissingSymbolicToolbox', ...
        'This verification script requires Symbolic Math Toolbox.');
end
if nargin<1 || isempty(path)
    path=fullfile(fileparts(mfilename('fullpath')),'exclusion_certificate.json');
end
cert=read_certificate(path);
require(strcmp(cert.format,'rhomboid-forbidden-region-exclusion-v1'),'Unknown format.');
root=rational_tuple(cert.root);
require(same_tuple(root,{exact(137,512),exact(37,64),exact(137,512),exact(1)}),'Wrong root rectangle.');
a=root{1}; b=root{2}; c=root{3}; d=root{4};
require(0<a & a<b & b<1 & c==a & d==1,'Invalid root rectangle.');
require(a*a-4*a+1>0 & 3*b*b>1,'Root fails to enclose Omega.');
corner=rational_tuple(cert.corner); cornerSummary=verify_corner(corner);
paths=cell(1,numel(cert.leaves)); leafBoxes=cell(size(paths));
counts=containers.Map('KeyType','char','ValueType','double');
margins=containers.Map('KeyType','char','ValueType','any');
for index=1:numel(cert.leaves)
    leaf=cert.leaves(index); box=rational_tuple(leaf.box);
    require(numel(box)==4,'Invalid leaf box size.');
    require([box{1}<box{2};box{3}<box{4}],'Invalid leaf box.');
    paths{index}=leaf.path; leafBoxes{index}=box;
    ul=box{1}; uh=box{2}; vl=box{3}; vh=box{4}; reason=leaf.reason;
    if strcmp(reason,'outside_v_le_u')
        require(vh<=ul,'Incorrect v<=u exclusion.');
    elseif any(strcmp(reason,{'outside_u_low','outside_u_high','outside_forbidden'}))
        p=geometric_predicates(iv(ul,uh),iv(vl,vh));
        if strcmp(reason,'outside_forbidden')
            require(p.DEI.hi<0 & p.Hregion.hi<0 & (p.Ctheta.hi<0 | p.Cdist.hi<0), ...
                'Incorrect forbidden-region exclusion.');
        else
            predicate=p.(reason);
            require(predicate.lo>=0,'Incorrect u-domain exclusion.');
        end
    elseif any(strcmp(reason,{'H+','H-'}))
        value=evaluate_leaf(box,'H');
        if strcmp(reason,'H+'), margin=value.lo; else, margin=-value.hi; end
        require(margin>0,sprintf('H sign failed at leaf %d.',index-1));
        if isKey(margins,reason), margin=min([margins(reason),margin]); end
        margins(reason)=margin;
    elseif strcmp(reason,'corner')
        require(corner{1}<=ul & ul<uh & uh<=corner{2} & ...
            corner{3}<=vl & vl<vh & vh<=corner{4}, ...
            'Corner leaf lies outside certified rectangle.');
    else
        error('kite:UnknownReason','Unknown leaf reason: %s',reason);
    end
    if isKey(counts,reason), counts(reason)=counts(reason)+1; else, counts(reason)=1; end
end
count=verify_partition(paths,leafBoxes,root);
fprintf('PASS: %d exact rational rectangles form a complete partition.\n',count);
print_certificate_counts(counts,'Counts');
reasons=sort(keys(margins));
for k=1:numel(reasons)
    fprintf('%s: every signed endpoint exceeds %s\n', ...
        reasons{k},decimal(margins(reasons{k}),false,12));
end
fprintf('THEOREM VERIFIED: H=0 misses all points of Omega satisfying\n');
fprintf('  v*(1-u^2)-2*u >= 0, OR u^2+2*u*v-1 >= 0, OR\n');
fprintf('  [u^2+2*u-1 >= 0 AND 1-2*u*v-v^2 >= 0].\n');
fprintf('Therefore Z intersects none of the relative closures C,D,E,H,I.\n');
fprintf('No mass-ordering or mass-positivity theorem was used in this certificate.\n');
summary=certificate_summary(count,counts,margins); summary.corner=cornerSummary;
end

function summary = verify_corner(corner)
require(same_tuple(corner,{exact(57,100),exact(29,50),exact(57,100),exact(3,5)}),'Wrong corner box.');
ul=corner{1}; uh=corner{2}; vl=corner{3}; vh=corner{4};
u=iv(ul,uh); v=iv(vl,vh);
require(ul>exact(1,2) & ul==vl,'Corner comparison inequalities failed.');
require(3*ul*ul<1 & 1<3*uh*uh,'Corner must enclose u=1/sqrt(3).');
p=geometric_predicates(u,v);
require(p.DEI.hi<0,'DEI not excluded in corner.');
hvBounds=cell(1,4); index=0;
for i=0:1
    for j=0:1
        subbox={ul+(uh-ul)*exact(i,2),ul+(uh-ul)*exact(i+1,2), ...
            vl+(vh-vl)*exact(j,2),vl+(vh-vl)*exact(j+1,2)};
        value=evaluate_leaf(subbox,'Hv');
        require(value.hi<-1,'H_v<-1 failed in corner.');
        index=index+1; hvBounds{index}=value;
    end
end
lo=hvBounds{1}.lo; hi=hvBounds{1}.hi;
for i=2:4, lo=min([lo,hvBounds{i}.lo]); hi=max([hi,hvBounds{i}.hi]); end
hvHull=iv(lo,hi);
hmax=1-3*ul*ul;
require(uh+exact(1,2)*hmax<vh,'C comparison path leaves the corner box.');
vH=divide(subtract(1,multiply(u,u)),multiply(2,u));
require(vl<vH.lo & vH.lo<=vH.hi & vH.hi<vh,'H boundary path leaves the corner box.');
qc=c_path_quotient(u);
require(qc.lo>exact(3,5),'Positive C-path quotient failed.');
[b,K,E]=h_boundary_factors(u);
require(b.lo>0 & K.lo>0,'H-boundary positive factors failed.');
require(E.hi<-exact(3,8),'Negative H-boundary factor failed.');
fprintf('PASS corner: H_v subset %s and H_v<-1.\n',show(hvHull));
fprintf('PASS corner: Q_C subset %s and Q_C>3/5.\n',show(qc));
fprintf('PASS corner: E subset %s and E<-3/8.\n',show(E));
fprintf('PASS corner: both comparison paths remain in the corner rectangle.\n');
summary=struct('Hv',show(hvHull),'Q_C',show(qc),'E',show(E));
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

function c = divide(a,b)
c=multiply(a,powerExact(b,-1));
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

function yes = same_tuple(a,b)
yes=numel(a)==numel(b);
if ~yes, return; end
for k=1:numel(a)
    if ~proved(a{k}==b{k}), yes=false; return; end
end
end

function cert = read_certificate(path)
raw=fileread(path);
unquoted=regexprep(raw,'"(?:[^"\\]|\\.)*"','""');
require(isempty(regexp(unquoted, ...
    '(^|[\s\[\{,:])-?[0-9]+(\.[0-9]+([eE][+-]?[0-9]+)?|[eE][+-]?[0-9]+)', ...
    'once')), 'Floating-point number in exact certificate.');
cert=jsondecode(raw);
require(isstruct(cert) && isscalar(cert) && isfield(cert,'format') && ...
    isfield(cert,'root') && isfield(cert,'leaves'),'Malformed certificate.');
require(isstruct(cert.leaves) && ~isempty(cert.leaves),'Missing certificate leaves.');
end

function count = verify_partition(paths,leafBoxes,expectedRoot)
% Reconstruct the entire binary partition using exact symbolic comparisons.
require(iscell(paths) && iscell(leafBoxes) && numel(paths)==numel(leafBoxes) && ...
    ~isempty(paths),'Invalid partition data.');
for k=1:numel(paths)
    word=paths{k};
    require(ischar(word) && (isrow(word) || isempty(word)) && ...
        all(word=='0' | word=='1'),'Invalid binary path.');
end
capacity=1+sum(cellfun(@numel,paths));
children=zeros(capacity,2); boxes=cell(capacity,1); used=1;
for k=1:numel(paths)
    box=leafBoxes{k};
    require(iscell(box) && numel(box)==4,'Invalid leaf rectangle.');
    require([box{1}<box{2};box{3}<box{4}],'Degenerate or reversed rectangle.');
    node=1;
    for digit=paths{k}
        require(isempty(boxes{node}),'A leaf is an ancestor of another leaf.');
        branch=double(digit-'0')+1; child=children(node,branch);
        if child==0
            used=used+1; child=used; children(node,branch)=child;
        end
        node=child;
    end
    require(isempty(boxes{node}) && all(children(node,:)==0), ...
        'Duplicate leaf or conflicting path prefix.');
    boxes{node}=box;
end
counts=zeros(used,1);
for node=used:-1:1
    pair=children(node,:);
    if ~isempty(boxes{node})
        require(all(pair==0),'Leaf also has children.'); counts(node)=1;
        continue
    end
    require(all(pair>0),'Incomplete binary partition.');
    left=boxes{pair(1)}; right=boxes{pair(2)};
    vertical=proved(left{2}==right{1}) && same_tuple(left(3:4),right(3:4));
    horizontal=proved(left{4}==right{3}) && same_tuple(left(1:2),right(1:2));
    require(vertical || horizontal,'Children fail to partition their parent.');
    if vertical
        boxes{node}={left{1},right{2},left{3},left{4}};
    else
        boxes{node}={left{1},left{2},left{3},right{4}};
    end
    counts(node)=counts(pair(1))+counts(pair(2));
end
count=counts(1);
require(same_tuple(boxes{1},expectedRoot), ...
    'Partition does not cover the root rectangle exactly.');
require(count==numel(paths),'Partition leaf count mismatch.');
end

function value = evaluate_leaf(box,quantity)
u=iv(box{1},box{2}); v=iv(box{3},box{4});
if strcmp(quantity,'H'), value=mass_difference_numerator(u,v); return; end
require(any(strcmp(quantity,{'Hu','Hv'})),'Unknown derivative.');
result=mass_difference_numerator(dv(u,double(strcmp(quantity,'Hu'))), ...
    dv(v,double(strcmp(quantity,'Hv'))));
value=result.derivative;
end

function summary = certificate_summary(count,counts,margins)
reasons=sort(keys(counts)); rows=struct('reason',{},'count',{});
for k=1:numel(reasons), rows(k)=struct('reason',reasons{k},'count',counts(reasons{k})); end
names=sort(keys(margins)); bounds=struct('reason',{},'lower_bound',{});
for k=1:numel(names)
    bounds(k)=struct('reason',names{k},'lower_bound',decimal(margins(names{k}),false,12));
end
summary=struct('leaf_count',count,'counts',rows,'signed_endpoint_lower_bounds',bounds);
end

function print_certificate_counts(counts,label)
reasons=sort(keys(counts)); parts=cell(1,numel(reasons));
for k=1:numel(reasons), parts{k}=sprintf('"%s": %d',reasons{k},counts(reasons{k})); end
fprintf('%s: {%s}\n',label,strjoin(parts,', '));
end

function result = geometric_predicates(u,v)
% Exact interval/dual arithmetic in the original evaluation order.
result=struct('outside_u_low',add(subtract(multiply(u,u),multiply(4,u)),1),'outside_u_high', ...
    subtract(multiply(multiply(3,u),u),1),'DEI',subtract(multiply(v,subtract(1,multiply(u,u))), ...
    multiply(2,u)),'Hregion',subtract(add(multiply(u,u),multiply(multiply(2,u),v)),1),'Ctheta', ...
    subtract(add(multiply(u,u),multiply(2,u)),1),'Cdist',subtract(subtract(1,multiply(multiply(2,u), ...
    v)),multiply(v,v)));
end

function result = c_path_quotient(u)
% Exact interval/dual arithmetic in the original evaluation order.
h=subtract(1,multiply(multiply(3,u),u));
k=exact(1,2);
v=add(u,multiply(k,h));
K=divide(add(subtract(multiply(3,powerExact(u,4)),multiply(multiply(6,u),u)),7),powerExact(add(1, ...
    multiply(u,u)),3));
av=divide(v,add(1,multiply(v,v)));
au=divide(u,add(1,multiply(u,u)));
KC=multiply(divide(multiply(divide(multiply(8,powerExact(subtract(1,multiply(u,u)),3)),powerExact(u, ...
    3)),subtract(1,multiply(u,v))),multiply(add(1,multiply(v,v)),add(1,multiply(u,u)))), ...
    add(add(multiply(av,av),multiply(av,au)),multiply(au,au)));
T=divide(multiply(k,add(1,multiply(u,v))),multiply(multiply(2,v),subtract(1,multiply(u,u))));
t=multiply(h,T);
b=subtract(1,t);
q=multiply(t,subtract(1,t));
A=divide(multiply(8,powerExact(subtract(1,multiply(u,u)),3)),powerExact(add(1,multiply(u,u)),3));
B=divide(powerExact(subtract(1,multiply(u,u)),3),multiply(8,powerExact(u,3)));
C=divide(multiply(multiply(8,powerExact(v,3)),powerExact(subtract(1,multiply(u,u)),3)), ...
    multiply(powerExact(u,3),powerExact(add(1,multiply(v,v)),3)));
N=add(subtract(A,B),multiply(multiply(q,q),subtract(subtract(multiply(A,A),B),multiply(subtract(A,1),C))));
LT=add(subtract(multiply(multiply(b,b),K),multiply(T,add(add(multiply(b,b),b),1))), ...
    multiply(powerExact(b,3),add(K,multiply(k,KC))));
result=subtract(multiply(LT,N),multiply(multiply(multiply(multiply(multiply(b,b),subtract(A,B)),k), ...
    KC),subtract(1,q)));
end

function [b,K,E] = h_boundary_factors(u)
% Exact interval/dual arithmetic in the original evaluation order.
h=subtract(1,multiply(multiply(3,u),u));
K=divide(add(subtract(multiply(3,powerExact(u,4)),multiply(multiply(6,u),u)),7),powerExact(add(1, ...
    multiply(u,u)),3));
A=divide(multiply(8,powerExact(subtract(1,multiply(u,u)),3)),powerExact(add(1,multiply(u,u)),3));
B=divide(powerExact(subtract(1,multiply(u,u)),3),multiply(8,powerExact(u,3)));
b=divide(powerExact(add(1,multiply(u,u)),2),multiply(4,powerExact(subtract(1,multiply(u,u)),2)));
Q0=divide(multiply(subtract(3,multiply(u,u)),powerExact(add(1,multiply(u,u)),2)),multiply(16, ...
    powerExact(subtract(1,multiply(u,u)),4)));
E=add(multiply(subtract(A,B),add(negate(K),multiply(A,Q0))),multiply(multiply(h,powerExact(Q0,2)), ...
    subtract(subtract(multiply(multiply(2,A),A),powerExact(A,3)),B)));
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

function value = dv(a,derivative)
% Interval value and enclosure of one partial derivative.
if isDual(a), value=a; return; end
if nargin<2, derivative=0; end
value=struct('value',iv(a),'derivative',iv(derivative));
end

function yes = isDual(a)
yes=isstruct(a) && isfield(a,'value') && isfield(a,'derivative');
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

function result = mass_difference_numerator(u,v)
% Exact interval/dual arithmetic in the original evaluation order.
t=divide(multiply(subtract(v,u),add(1,multiply(u,v))),multiply(multiply(2,v),subtract(1,multiply(u,u))));
b=subtract(1,t);
q=multiply(t,subtract(1,t));
A=divide(multiply(8,powerExact(subtract(1,multiply(u,u)),3)),powerExact(add(1,multiply(u,u)),3));
B=divide(powerExact(subtract(1,multiply(u,u)),3),multiply(8,powerExact(u,3)));
C=divide(multiply(multiply(8,powerExact(v,3)),powerExact(subtract(1,multiply(u,u)),3)), ...
    multiply(powerExact(u,3),powerExact(add(1,multiply(v,v)),3)));
N=add(subtract(A,B),multiply(multiply(q,q),subtract(subtract(multiply(A,A),B),multiply(subtract(A,1),C))));
L=subtract(add(multiply(multiply(b,b),subtract(A,1)),multiply(powerExact(b,3),C)),1);
result=subtract(multiply(L,N),multiply(multiply(multiply(multiply(b,b),subtract(A,B)),subtract(C, ...
    A)),subtract(1,q)));
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
