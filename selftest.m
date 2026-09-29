function selftest(varargin)
%SELFTEST Test symbolic rationals, interval bounds and certificate rejection.
%   selftest                 % fast tests
%   selftest('Full',true)     % all certificates + six reference comparisons
if isempty(ver('symbolic'))
    error('selftest:MissingSymbolicToolbox', ...
        'This verification script requires Symbolic Math Toolbox.');
end
p=inputParser;
addParameter(p,'Full',false,@(v)islogical(v)&&isscalar(v));
parse(p,varargin{:}); base=fileparts(mfilename('fullpath'));
require(exact('0.31')==sym(31)/100,'Exact decimal parsing');
require(exact('9007199254740993')-exact('9007199254740992')==1,'Large integer parsing');
require(exact('6/-8')==sym(-3)/4,'Rational normalization');
require(exact(intmax('uint64'))==sym('18446744073709551615'),'Integer input precision');
require(strcmp(decimal(sym(-1)/3,false,4),'-0.3334') && ...
    strcmp(decimal(sym(-1)/3,true,4),'-0.3333'),'Negative directed rounding');
require(strcmp(decimal(sym(-1)/100000,true,4),'0.0000'),'Directed rounding across zero');
require(strcmp(decimal(sym(1)/3,false,0),'0') && strcmp(decimal(sym(1)/3,true,0),'1'), ...
    'Integer directed rounding');
cases=jsondecode(fileread(fullfile(base,'arithmetic_reference.json')));
for k=1:numel(cases)
    c=cases(k); a=exact(c.a); b=exact(c.b);
    equalZero([a+b-exact(c.add);a-b-exact(c.subtract);a*b-exact(c.multiply); ...
        a/b-exact(c.divide);a^c.power-exact(c.powered)], ...
        sprintf('Python Fraction arithmetic case %d',k));
    require(sign(a-b)==c.comparison,'Rational comparison');
    require(strcmp(decimal(a,false,12),c.lower) && strcmp(decimal(a,true,12),c.upper), ...
        'Python Fraction decimal reference');
end
syms t real
equalZero((t+1)^2-(t^2+2*t+1),'Exact polynomial identity');
expect_error(@()require(t>0,'Undecidable statement'),'kite:CheckFailed');
i=iv(-2,3); squared=powerExact(i,2); odd=powerExact(i,3);
require([squared.lo==0;squared.hi==9;odd.lo==-8;odd.hi==27],'Interval powers');
negative=powerExact(iv(-4,-2),-1);
require([negative.lo==sym(-1)/2;negative.hi==sym(-1)/4],'Negative interval reciprocal');
product=multiply(iv(-2,3),iv(-4,5));
require([product.lo==-12;product.hi==15],'Four-endpoint interval product');
mixed={add(sym(1)/2,iv(1,2)),add(iv(1,2),sym(1)/2)};
for k=1:2, require([mixed{k}.lo==sym(3)/2;mixed{k}.hi==sym(5)/2],'Scalar/interval addition'); end
d=dv(iv(2),1);
mixed={multiply(sym(1)/2,d),multiply(d,sym(1)/2), ...
    multiply(iv(sym(1)/2),d),multiply(d,iv(sym(1)/2))};
for k=1:4
    require(isDual(mixed{k}),'Mixed derivative type');
    require([mixed{k}.value.lo==1;mixed{k}.value.hi==1; ...
        mixed{k}.derivative.lo==sym(1)/2;mixed{k}.derivative.hi==sym(1)/2], ...
        'Mixed derivative arithmetic');
end
f=powerExact(d,3);
require([f.value.lo==8;f.value.hi==8;f.derivative.lo==12;f.derivative.hi==12],'Dual power');
f=divide(1,d);
require([f.value.lo==sym(1)/2;f.value.hi==sym(1)/2; ...
    f.derivative.lo==sym(-1)/4;f.derivative.hi==sym(-1)/4],'Dual reciprocal');
expect_error(@()exact(0.1),'kite:CheckFailed');
expect_error(@()exact(1,0),'kite:CheckFailed');
expect_error(@()iv(2,1),'kite:CheckFailed');
expect_error(@()divide(1,iv(-1,1)),'kite:CheckFailed');
root={sym(0),sym(2),sym(0),sym(1)};
left={sym(0),sym(1),sym(0),sym(1)}; right={sym(1),sym(2),sym(0),sym(1)};
require(verify_partition({'0','1'},{left,right},root)==2,'Valid binary partition');
expect_error(@()verify_partition({'0','0'},{left,left},root),'kite:CheckFailed');
expect_error(@()verify_partition({'0','00'},{left,left},root),'kite:CheckFailed');
expect_error(@()verify_partition({'0',''},{left,root},root),'kite:CheckFailed');
expect_error(@()verify_partition({'x'},{root},root),'kite:CheckFailed');
expect_error(@()verify_partition({'0'},{left},root),'kite:CheckFailed');
gap={sym(11)/10,sym(2),sym(0),sym(1)};
expect_error(@()verify_partition({'0','1'},{left,gap},root),'kite:CheckFailed');
wrongRoot={sym(0),sym(3),sym(0),sym(1)};
expect_error(@()verify_partition({'0','1'},{left,right},wrongRoot),'kite:CheckFailed');
temporary=[tempname '.json']; cleanup=onCleanup(@()delete_if_present(temporary));
write_utf8(temporary,'{"format":"test","root":[0.31],"leaves":[]}');
expect_error(@()read_certificate(temporary),'kite:CheckFailed');
write_utf8(temporary,'{"format":"test","root":[1e3],"leaves":[]}');
expect_error(@()read_certificate(temporary),'kite:CheckFailed');
fprintf('PASS: symbolic arithmetic, outward rounding, derivatives and invalid-certificate tests.\n');
if p.Results.Full
    report=run_all(); expected=jsondecode(fileread(fullfile(base,'reference_expected.json')));
    expectedIds=fieldnames(expected); actualIds={report.checks.id};
    require(all(ismember(expectedIds,actualIds)),'A reference certificate was not run.');
    for k=1:numel(expectedIds)
        id=expectedIds{k}; index=find(strcmp(actualIds,id),1);
        actual=jsondecode(jsonencode(report.checks(index).details));
        require(same_json(actual,expected.(id)),['Python reference mismatch: ' id]);
    end
    fprintf('PASS: all %d certificates completed; %d summaries equal the Python reference outputs.\n', ...
        numel(report.checks),numel(expectedIds));
end
end

function expect_error(f,expectedId)
try
    f();
catch exception
    require(strcmp(exception.identifier,expectedId), ...
        ['Unexpected failure: ' exception.identifier ' (expected ' expectedId ')']);
    return
end
error('kite:Selftest','An invalid operation was accepted.');
end

function yes = same_json(a,b)
if isstruct(a) && isstruct(b)
    names=sort(fieldnames(a)); yes=isequal(names,sort(fieldnames(b))) && numel(a)==numel(b);
    if ~yes, return; end
    for k=1:numel(a)
        for j=1:numel(names)
            if ~same_json(a(k).(names{j}),b(k).(names{j})), yes=false; return; end
        end
    end
elseif iscell(a) && iscell(b)
    yes=numel(a)==numel(b);
    if ~yes, return; end
    for k=1:numel(a), if ~same_json(a{k},b{k}), yes=false; return; end, end
elseif isnumeric(a) && isnumeric(b)
    yes=isequal(a(:),b(:));
else
    yes=isequal(a,b);
end
end

function delete_if_present(path)
if exist(path,'file')==2, delete(path); end
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

function write_utf8(path,text)
fid=fopen(path,'w','n','UTF-8');
require(fid>=0,['Cannot write file: ' path]);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',text);
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

function yes = same_tuple(a,b)
yes=numel(a)==numel(b);
if ~yes, return; end
for k=1:numel(a)
    if ~proved(a{k}==b{k}), yes=false; return; end
end
end
