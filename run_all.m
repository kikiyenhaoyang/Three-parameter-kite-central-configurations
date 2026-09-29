function report = run_all(varargin)
%RUN_ALL Verify all bundled exact certificates and save logs plus summary.json.
%   run_all
%   run_all('Only',{'K0'})
%   run_all('Only',{'S1','B1'})
%   run_all('Only',{'K2','R3'})
%   run_all('List',true)
%   run_all('CheckOnly',true)
%   report = run_all('OutputDir','my_results')
%
% Requires Symbolic Math Toolbox. All helpers are local to this file.
% Mathematical certificates use exact symbolic rationals; decimals are display only.
% Checks run as separate function calls in the current MATLAB process.
p=inputParser;
addParameter(p,'Only',{},@(v)iscellstr(v)||ischar(v)||isstring(v));
addParameter(p,'List',false,@(v)islogical(v)&&isscalar(v));
addParameter(p,'CheckOnly',false,@(v)islogical(v)&&isscalar(v));
addParameter(p,'OutputDir','',@(v)ischar(v)||(isstring(v)&&isscalar(v)));
parse(p,varargin{:}); args=p.Results;
if isempty(ver('symbolic'))
    error('run_all:MissingSymbolicToolbox', ...
        'This verification script requires Symbolic Math Toolbox.');
end
require(usejava('jvm'),'File checksum checking requires the MATLAB JVM.');
base=fileparts(mfilename('fullpath'));
files=check_bundle_files(base);
manifest=jsondecode(fileread(fullfile(base,'manifest.json')));
require(strcmp(manifest.format,'kite-matlab-symbolic-certificates-v2'),'Unknown manifest format.');
checks=manifest.checks; ids={checks.id};
require(numel(unique(ids))==numel(ids),'Duplicate certificate ID.');
for k=1:numel(checks)
    data=checks(k).data;
    if isstring(data), data=cellstr(data); end
    if isempty(data), data={}; end
    needed=[{checks(k).script};data(:)];
    for j=1:numel(needed)
        require(any(strcmp(needed{j},files)), ...
            ['A certificate file is absent from SHA256SUMS: ' needed{j}]);
    end
    [~,functionName,ext]=fileparts(checks(k).script);
    require(strcmp(ext,'.m') && strcmp(checks(k).entrypoint,functionName), ...
        'Manifest function/script mismatch.');
    expected=bundle_path(base,checks(k).script);
    actual=which(functionName);
    require(strcmp(actual,expected),['Another MATLAB function shadows ' functionName]);
end
report=struct('bundle_version',manifest.version,'checks',[]);
if args.List
    for k=1:numel(checks), fprintf('%s: %s\n',checks(k).id,checks(k).title); end
    return
end
if args.CheckOnly
    fprintf('PASS: certificate scripts, data and MATLAB support files are present.\n');
    return
end
only=args.Only;
if ischar(only), only={only}; elseif isstring(only), only=cellstr(only); end
if isempty(only)
    selected=checks;
else
    unknown=setdiff(only,ids);
    require(isempty(unknown),['Unknown certificate IDs: ' strjoin(unknown,', ')]);
    selected=checks(ismember(ids,only));
end
output=char(args.OutputDir);
if isempty(output), output=fullfile(base,'results'); end
outputFile=javaObject('java.io.File',output); output=char(outputFile.getCanonicalPath());
if exist(output,'dir')~=7
    [ok,message]=mkdir(output); require(ok,message);
end
report.started_utc=char(datetime('now','TimeZone','UTC', ...
    'Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
report.matlab=version; report.platform=computer;
report.arithmetic='Exact Symbolic Math Toolbox rational endpoints';
report.selected={selected.id};
report.checks=struct('id',{},'status',{},'returncode',{}, ...
    'elapsed_seconds',{},'log',{},'details',{});
totalStart=tic;
for k=1:numel(selected)
    item=selected(k);
    fprintf('RUN %s: %s\n',item.id,item.title); start=tic;
    [transcript,details,problem]=evalc('invoke_check(item.entrypoint)');
    elapsed=toc(start); passed=isempty(problem);
    if passed, status='PASS'; else, status='FAIL'; end
    logName=[item.id '.log'];
    write_utf8(fullfile(output,logName),transcript);
    report.checks(k)=struct('id',item.id,'status',status,'returncode',double(~passed), ...
        'elapsed_seconds',elapsed,'log',logName,'details',details);
    fprintf('%s %s (%.3f s)\n',status,item.id,elapsed);
end
report.elapsed_seconds=toc(totalStart);
report.all_passed=all(strcmp({report.checks.status},'PASS'));
write_utf8(fullfile(output,'summary.json'),[jsonencode(report) newline]);
if report.all_passed
    fprintf('ALL %d SELECTED CERTIFICATES PASSED.\n',numel(selected));
else
    error('kite:CertificateFailed','At least one certificate failed; see logs in %s',output);
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

function path = bundle_path(base,relative)
require(ischar(relative) && ~isempty(relative) && ...
    isempty(regexp(relative,'[^A-Za-z0-9_./-]','once')) && relative(1)~='/' && ...
    ~any(strcmp(strsplit(relative,'/'),'..')), 'Invalid bundle-relative path.');
baseFile=javaObject('java.io.File',base); canonicalBase=char(baseFile.getCanonicalPath());
file=javaObject('java.io.File',fullfile(base,relative)); path=char(file.getCanonicalPath());
prefix=[canonicalBase filesep];
if ispc, contained=strncmpi(path,prefix,numel(prefix));
else, contained=strncmp(path,prefix,numel(prefix)); end
require(contained,['A bundle path points outside the archive: ' relative]);
end

function write_utf8(path,text)
fid=fopen(path,'w','n','UTF-8');
require(fid>=0,['Cannot write file: ' path]);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',text);
end

function [details,problem] = invoke_check(functionName)
% Catch inside evalc so output preceding an error is retained in the log.
details=[]; problem=[];
try
    if nargout(functionName)==0
        feval(functionName);
    else
        details=feval(functionName);
    end
catch exception
    problem=exception;
    fprintf('\n[error]\n%s\n',getReport(exception,'extended','hyperlinks','off'));
end
end

function files = check_bundle_files(base)
lines=regexp(fileread(fullfile(base,'SHA256SUMS')),'\r?\n','split');
if ~isempty(lines) && isempty(lines{end}), lines(end)=[]; end
files=cell(1,numel(lines));
for k=1:numel(lines)
    tokens=regexp(lines{k},'^([0-9a-f]{64})  (.+)$','tokens','once');
    require(~isempty(tokens),'Malformed checksum entry.');
    relative=tokens{2};
    require(~any(strcmp(relative,files(1:k-1))),'Duplicate checksum entry.');
    path=bundle_path(base,relative);
    require(exist(path,'file')==2,['Missing bundle file: ' relative]);
    require(strcmp(sha256_file(path),tokens{1}),['Checksum mismatch: ' relative]);
    files{k}=relative;
end
required={'manifest.json','run_all.m','README.md','selftest.m','reference_expected.json'};
require(all(ismember(required,files)),'Incomplete checksum inventory.');
fprintf('PASS: SHA-256 checksums for %d files.\n',numel(files));
end

function digest = sha256_file(path)
fid=fopen(path,'rb');
require(fid>=0,['Cannot open file: ' path]);
cleanup=onCleanup(@()fclose(fid));
bytes=fread(fid,Inf,'*uint8');
md=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
md.update(typecast(bytes(:),'int8'));
raw=typecast(md.digest(),'uint8');
digest=lower(reshape(dec2hex(raw,2).',1,[]));
end
