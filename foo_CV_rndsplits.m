%
% foo_CV_rndsplits( 'Fundus_groups.dt','../_computed_sigs_29img' );
%
% By analogy with 'foo_Gleason_rndsplits.m', generate CV splits for mat files from
% the given path (mat_files_path) and given model (class information, eg ground truth). 
% The model is given in separate model file (Fundus_groups.dt).
%

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function foo_CV_rndsplits( varargin )
tic;
modelfile = varargin{1};
mats_path = varargin{2};

testingTYPE = 'image';

% A subject to reconstruct...
%[categSet,categ_CoreList] = read_model( modelfile );

if ~exist( mats_path,'dir' ), error(['**path ' mats_path 'doesn''t exist,abort**']);end 

% mat files here:
exten = 'mat';

foo = dir( fullfile( mats_path,['*.' exten] ) );
pool = {foo.name};
[unique_corenames,nhits,unique_core_letter] = get_tileCore( pool,exten );

% split mats into classes
categorized_pool = categorize_pool( pool,categSet,categ_CoreList,exten );
clear pool;

% make splits by randomizing classes
perclassTR = 400;
nCV = 30;
%perclassTR = 20;nCV = 2;

%perclassTR = 500;
%nCV = 2;

name_templ = modelfile;

switch testingTYPE
 case 'core'
  gen_perCORE_splits( unique_corenames,categorized_pool,nCV,perclassTR,name_templ,categSet,mats_path,exten );
 case 'image'
  gen_perimage_splits( categorized_pool,nCV,perclassTR,name_templ,categSet,mats_path );
    otherwise,fprintf('\n>>warning: testingTYPE given bans making splits<<\n');
end % switch

fprintf(':: %u splits are now complete ::\n',nCV);
toc;
end % eofunc


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function gen_perCORE_splits( unique_corenames,categorized_pool,nCV,NperclassTR,name_templ,categSet,mats_path,exten )

[foo,name_templ,ext] = fileparts( name_templ );

wrk = 'work';
mymkdir(wrk);
name_templ = fullfile(wrk,name_templ);

TRtempl = [name_templ '_TR'];
TEtempl = [name_templ '_TE'];
Ncla = length(categorized_pool);

for cla = 1:Ncla,
  uCORES_present{cla} = get_allperclass_Cores( categorized_pool{cla},exten );
end % cla

for sp = 1:nCV,
 fprintf('Split %2u/%u: ',sp,nCV);
 TRtempl_sp = [TRtempl sprintf('_Split%u_perCOR.mat',sp)];
 TEtempl_sp = [TEtempl sprintf('_Split%u_perCOR.mat',sp)];
 
 TR_paths = {};     TE_paths = {};
 image_ids_TR = []; image_ids_TE = [];
 class_vecTR = [];  class_vecTE = [];
 
 % random shaffling is performed for each category
 for cla = 1:Ncla,
  uniqueCORESperclass = uCORES_present{cla};
  Nuniq  = length(uniqueCORESperclass);
  rp_ind = randperm( Nuniq );
  uniqueCORESperclass = uniqueCORESperclass(rp_ind);
  
  [TR,TE,TR_imIDs,TE_imIDs] = comile_mats_PCORE_TRTEsplit( categorized_pool{cla},uniqueCORESperclass,NperclassTR );
  
  class_vecTR = [class_vecTR cla.*ones(1,length(TR))];
  class_vecTE = [class_vecTE cla.*ones(1,length(TE))];
  
  TR_paths = [TR_paths TR];
  TE_paths = [TE_paths TE];
  image_ids_TR = [image_ids_TR TR_imIDs];
  image_ids_TE = [image_ids_TE TE_imIDs];
  clear TR TE;
  
 end % cla
 
 % assemble both TR/TE sig matrices
 [TR,TE,signature_labels] = assemble_sigmatr(mats_path,TR_paths,TE_paths,class_vecTR,class_vecTE);
  
 category_names = categSet;
 dataset_name = name_templ;
 image_ids = image_ids_TR;
 image_paths = TR_paths;
 signature_matrix = TR;
 nmlist = {'category_names','dataset_name','image_paths','signature_labels','signature_matrix','image_ids'};
 mysave(TRtempl_sp,nmlist,category_names,dataset_name,image_paths,signature_labels,signature_matrix,image_ids);
 %%save(TRtempl_sp,'category_names','dataset_name','image_paths','signature_labels','signature_matrix','image_ids');

 image_ids = image_ids_TE;
 image_paths = TE_paths;
 signature_matrix = TE;
 classvecTE = TE(end,:);
 moveTEmats( TRtempl_sp,mats_path,TE_paths,image_ids_TE,category_names,dataset_name,signature_labels,classvecTE );
 
 fprintf('%s/%s saved\n',TRtempl_sp,TEtempl_sp);
end % split
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function moveTEmats( target_path,src_path,indiv_TE_paths,image_ids_TE,category_names,dataset_name,signature_labels,classvecTE )
for TEind = 1:length(indiv_TE_paths),
 loadpath = fullfile( src_path,indiv_TE_paths{TEind} );
 load(loadpath);
 image_ids = image_ids_TE(TEind);
 if TEind==1,
  [fold,name,ext] = fileparts( target_path );
  name = fullfile(fold,name);
  mymkdir( name );
 end % if TEind==1
 cpName = fullfile( name,indiv_TE_paths{TEind} );
 
 signature_matrix = [dummyDDB;classvecTE(TEind)];
 image_paths = dummyNm;
%%save( cpName,'signature_matrix','image_paths','signature_labels','image_ids','category_names','dataset_name' );
 %mysave( cpName,'signature_matrix','image_paths','signature_labels','image_ids','category_names','dataset_name' );
 
 nmlist = {'category_names','dataset_name','image_paths','signature_labels','signature_matrix','image_ids'};
 mysave(cpName,nmlist,category_names,dataset_name,image_paths,signature_labels,signature_matrix,image_ids);
end % TEind
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function mymkdir( name )
if exist(name,'dir'),return;end
stat = mkdir( name );
if ~stat,error(['>>mymkdir<< **cannot make folder ' name ',abort**']);end
end % eofunc

% this is to deal with the fact the earlier matlab cannot save into
% specified folder
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function mysave( panm,nmlist,category_names,dataset_name,image_paths,signature_labels,signature_matrix,image_ids )
[fld,nm] = fileparts( panm );
homedir = pwd;
if ~isempty(fld),cd(fld);else,error('**empty ''fld'',abort**');end
save(nm,'category_names','dataset_name','image_paths','signature_labels','signature_matrix','image_ids');
cd(homedir);
return;

n = length(nmlist);
ss = [];
for ii = 1:n,
 cur_name = nmlist{ii};
 sfoo = sprintf('''%s''',nmlist{ii});
 ss = [ss sfoo ','];
end % ii
ss = ss(1:end-1);
foo = ['save(nm,' ss ');'];
eval( foo );

cd(homedir);
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function uniqueCORES = get_allperclass_Cores( category_pool,exten )
n = length(category_pool);
uniqueCORES = {}; allCORES = {};
for ii = 1:n,
 core = get_corename( category_pool{ii},exten,'_Grade' );
 allCORES{ii} = core;
end % ii
uniqueCORES = unique( allCORES );
end % eofunc

% Make per-image CV splits, keeping in each class number of samples (tiles)
% being about value 'perclassTR'. This number is not exact, because it
% could interfere with unfinished image: I finish the ongoing image, then
% break the class
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function gen_perimage_splits( categorized_pool,nCV,NperclassTR,name_templ,categSet,mats_path )

[foo,name_templ,ext] = fileparts( name_templ );

wrk = 'work';
mymkdir(wrk);
name_templ = fullfile(wrk,name_templ);

TRtempl = [name_templ '_TR'];
TEtempl = [name_templ '_TE'];
Ncla = length(categorized_pool);
for cla = 1:Ncla,
  clalen(cla) = length( categorized_pool{cla} );
  [uIDs{cla},imageIDs{cla}] = get_allperclass_imageIDs(categorized_pool{cla});
end % cla

for sp = 1:nCV,
 fprintf('Split %2u/%u: ',sp,nCV);
 TRtempl_sp = [TRtempl sprintf('_Split%u_perIm.mat',sp)];
 TEtempl_sp = [TEtempl sprintf('_Split%u_perIm.mat',sp)];
 
 TR_paths = {};     TE_paths = {};
 image_ids_TR = []; image_ids_TE = [];
 class_vecTR = [];  class_vecTE = [];
 
 % random shaffling is performed for each category
 for cla = 1:Ncla,
  uniqueIDs_perClass = uIDs{cla};
  all_imIDs_perClass = imageIDs{cla};
  
  Nuniq  = length(uniqueIDs_perClass);
  rp_ind = randperm(Nuniq);
  uniqueIDs_perClass = uniqueIDs_perClass(rp_ind);
  
  [TR,TE,TR_imIDs,TE_imIDs] = comile_mats_PI_TRTEsplit( categorized_pool{cla},uniqueIDs_perClass,NperclassTR );
  
  class_vecTR = [class_vecTR cla.*ones(1,length(TR))];
  class_vecTE = [class_vecTE cla.*ones(1,length(TE))];
  
  TR_paths = [TR_paths TR];
  TE_paths = [TE_paths TE];
  image_ids_TR = [image_ids_TR TR_imIDs];
  image_ids_TE = [image_ids_TE TE_imIDs];
  clear TR TE;
  
 end % cla
 
 % assemble both TR/TE sig matrices
 [TR,TE,signature_labels] = assemble_sigmatr(mats_path,TR_paths,TE_paths,class_vecTR,class_vecTE);
  
 category_names = categSet;
 dataset_name = name_templ;
 image_ids = image_ids_TR;
 image_paths = TR_paths;
 signature_matrix = TR;
%%save(TRtempl_sp,'category_names','dataset_name','image_paths','signature_labels','signature_matrix','image_ids');
 %mysave(TRtempl_sp,'category_names','dataset_name','image_paths','signature_labels','signature_matrix','image_ids');

 nmlist = {'category_names','dataset_name','image_paths','signature_labels','signature_matrix','image_ids'};
 mysave(TRtempl_sp,nmlist,category_names,dataset_name,image_paths,signature_labels,signature_matrix,image_ids);
 
 image_ids = image_ids_TE;
 image_paths = TE_paths;
 signature_matrix = TE;
 classvecTE = TE(end,:);
 moveTEmats( TRtempl_sp,mats_path,TE_paths,image_ids,category_names,dataset_name,signature_labels,classvecTE );
 
 fprintf('%s/%s saved\n',TRtempl_sp,TEtempl_sp);
end % split
end % eofunc


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function [TRmat,TEmat,siglabs] = assemble_sigmatr(mats_path,TR,TE,class_vecTR,class_vecTE)
TRmat = []; TEmat = [];
% TR matrix
for sa = 1:length(TR),
 curr = fullfile( mats_path,TR{sa} );
 load(curr);
%tmp1 = dummyDDB(end);
%tmp2 = class_vecTR(sa);
%if dummyDDB(end) ~= class_vecTR(sa),error('**mismatch category index,abort**');end
 TRmat = [TRmat dummyDDB];
end % sa
TRmat = [TRmat; class_vecTR];

% TE matrix
for sa = 1:length(TE),
 curr = fullfile( mats_path,TE{sa} );
 load(curr);
 TEmat = [TEmat dummyDDB];
end % sa
TEmat = [TEmat; class_vecTE];
end % eofunc

% compile sets of mat files, per-TR and per-TE based on PER-IMAGE 
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function [TR,TE,TR_imageIDs,TE_imageIDs] = comile_mats_PCORE_TRTEsplit( class_pool,queryCORES,NperclassTR )
TR = {}; TE = {};
TR_imageIDs = [];
TE_imageIDs = [];
% Fill TE with first CORE: get all samples for given core (that is always
% the very first core from the query)
TEcore_query = queryCORES{1};
matching_indic = get_COREmatches( TEcore_query,class_pool ); 
matching_names = class_pool(matching_indic);
TE = [TE matching_names];
TE_imageIDs = getIDs_for_matching_names( matching_names );
iLast = 1;

% Fill TR with the rest CORES from the query
for ii = iLast+1:length(queryCORES),
 TRcore_query = queryCORES{ii};
 matching_indic = get_COREmatches( TRcore_query,class_pool ); 
 matching_names = class_pool(matching_indic);
 if length(TR) >= NperclassTR,break;end
 if length(TR) + length(matching_names) > NperclassTR,
  safeNum = NperclassTR - length(TR);
  if safeNum < 0,error('** negative safeNum,abort **');end
  foo = matching_names(1:safeNum);
  TR = [TR foo]; 
  TR_imageIDs = [TR_imageIDs getIDs_for_matching_names( foo )];
  break;
 end % if
 TR = [TR matching_names];
 TR_imageIDs = [TR_imageIDs getIDs_for_matching_names( matching_names )];
end % ii
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function imageIDs = getIDs_for_matching_names( matching_names )
for ii = 1:length(matching_names),
 imageIDs(ii) = get_imID( matching_names{ii} );
end % ii
end % eofunc

% compile sets of mat files, per-TR and per-TE based on PER-IMAGE 
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function [TR,TE,TR_imageIDs,TE_imageIDs] = comile_mats_PI_TRTEsplit( class_pool,queryIDs,NperclassTR )
TR = {}; TE = {};
TR_imageIDs = [];
TE_imageIDs = [];
% Fill TR with first NperclassTR elements
for ii = 1:length(queryIDs),
 matching_indic = get_IDmatches ( queryIDs(ii),class_pool );
 matching_names = class_pool(matching_indic);
 TR = [TR matching_names];
 TR_imageIDs = [TR_imageIDs queryIDs(ii).*ones(1,length(matching_names))];
 if length(TR) >= NperclassTR, iLast = ii; break; end
end % ii

% Fill TE with the rest
for ii = iLast+1:length(queryIDs),
 matching_indic = get_IDmatches ( queryIDs(ii),class_pool );
 matching_names = class_pool(matching_indic);
 TE = [TE matching_names];
 TE_imageIDs = [TE_imageIDs queryIDs(ii).*ones(1,length(matching_names))];
end % ii
end % eofunc

function printcell(fname,arr,varargin),
if nargin > 2,k0=varargin{0};k1=varargin{1};
else,k0=1;k1=length(arr);end
fid = fopen(fname,'wt');
for jj = k0:k1,
 fprintf(fid,'%s\n',arr{jj});
end % jj
fclose(fid);
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function [uIDs,imageIDs] = get_allperclass_imageIDs( class_pool )
n = length(class_pool);
uIDs = []; imageIDs = [];

% just a test!
%matching_names = get_IDmatches ( 695,class_pool );

for ii = 1:n,
 currID = get_imID (class_pool{ii});
 imageIDs = [imageIDs currID];
 if ~ismember(currID,uIDs),uIDs = [uIDs currID];end
end % ii
end % eofunc
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function ID = get_imID ( nm ),
beg = '001_';en = '__';
k1 = strfind(nm,en);
k0 = strfind(nm,beg);
foo = (k0+length(beg)):(k1-1);
ID = str2num( nm(foo) );
end % eofunc

% Get all file names corresponding to given Core name
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function matching_indices = get_COREmatches ( core_name,class_pool ),
ind = strfind( class_pool,sprintf('%s_',core_name) );
% found this at: http://arnabocean.com/frontposts/2013-10-14-matlabstrfind/
%matching_indices = find( ~cellfun(@isempty,ind) );
matching_indices = find( ~cellfun(@isempty,ind) & cellfun(@(x) isequal(x,1),ind) ); % must be in the beginning of str

% May 16:
% It was giving a bug while seeking for 'e2_' and was returning it at every
% occursion of '.._Grade2_Gleason6_E_Auto_188_50x-001...' -- it was not
% Grade2 safe!

%%%matching_names = class_pool(matching_names);
end % eofunc
function yn = eq1(x)
yn = logical(0);
if x==1,yn=logical(1);end
end % eofunc

% Get all file names corresponding to given image ID
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function matching_names = get_IDmatches ( ID,class_pool ),
ID = num2str(ID);
ind = strfind( class_pool,sprintf('_%s__',ID) );
% found this at: http://arnabocean.com/frontposts/2013-10-14-matlabstrfind/
matching_names = find( ~cellfun(@isempty,ind) );
end % eofunc

% Split the whole pool into categories
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function categorized_pool = categorize_pool( pool,categSet,categ_CoreList,exten )
n = length(pool);
categorized_pool = cell(1,length(categ_CoreList));
for ii = 1:n,
 %fprintf('categorizing progress: %5u/%u\n',ii,n);
 curr = pool{ii};
 catg = assign_class( curr,categ_CoreList,exten );
 if isempty(catg),continue;end
 categorized_pool{catg}{end+1} = curr;
end % c
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function catg = assign_class( curr,categ_CoreList,exten )
core = get_corename( curr,exten,'_Grade' ); catg = [];
Ncla = length( categ_CoreList );
for c = 1:Ncla,
 if ismember( core,categ_CoreList{c} ),catg = c;break;end
end % c
% catg may be empty: no match means the core is out of the assigned set,
% this is allowed
if isempty(catg),return;end
if ~ismember(catg,1:Ncla),error('>>assign_class<< **out of range category,abort**');end
end % eofunc

% Read CategoryNames and Category core list
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function [categSet,categ_CoreList] = read_model( modelfile )
if ~exist(modelfile,'file'),error(['**model ' modelfile ' doesn''t exis,abort**']);end
fid = fopen( modelfile ); foo = fgetl(fid);
while ~strcmpi( foo,'#EOH' ),foo = fgetl(fid);end
cnt = 0;
while ~strcmpi( foo,'#EOF' ),
 foo = fgetl(fid);
 if strcmpi( foo,'#EOF' ),continue;end
 kname = strfind( foo,':' );if isempty(kname),error('**category name problem,abort**');end
 cnt = cnt+1;
 categSet{cnt} = foo(1:(kname(1)-1));
 k0 = strfind( foo,'{' );
 k1 = strfind( foo,'}' );
 if isempty(k0)|isempty(k1),error('**category core_list problem,abort**');end
 corelist = foo((k0(1)+1):(k1(1)-1));
 C = textscan( corelist,'%s','delimiter',',');
 categ_CoreList{cnt} = C{1}';
end % while
fclose(fid);
end % eofunc


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function [unique_corenames,nhits,unique_core_letter] = get_tileCore( pool,exten )

pattn = ['.' exten];
cnt = 0;
for ii = 1:length(pool),
 curr = pool{ii};
 core = get_corename( curr,exten,'_Grade' );
 cnt = cnt + 1;
 wholepool{cnt} = core;
end % ii

unique_corenames = unique( wholepool );
wholepool = sort( wholepool );

% compute number of hits per core
for c = 1:length(unique_corenames),u = unique_corenames{c};
 foo = find(strcmpi(u,wholepool));
 nhits(c) = length(foo);
 unique_core_letter{c} = u(1);
end % c
unique_core_letter = unique(unique_core_letter);

report_data(unique_corenames,nhits,unique_core_letter);
end % eofunc
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function core = get_corename( str,exten,suffix )
pattn = ['.' exten]; k = strfind( str,pattn );
if isempty(k),error(['>>get_corename<<' ' **pattn ' pattn ' not found in given name,abort**']);end
nm = str(1:(k(1)-1));
kCore = strfind( nm,'_Grade' ); 
if isempty(kCore),error(['>>get_corename<<' ' ** ' '_Grade' ' not found in given name,abort**']);end
core = nm(1:(kCore-1));
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function report_data(unique_corenames,nhits,unique_core_letter)
fprintf('\nReport cores and per-core hits\n');
fprintf('( eg number of tiles computed )\n\n');
for le = unique_core_letter,
 core_instances = mystrfind( unique_corenames,le{1} );
 buff = unique_corenames(core_instances);
 les  = nhits(core_instances);
 printem( buff,les );
end % le
fprintf('\n');
end % eofunc
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function printem( buff,les )
for ii = 1:length(buff),
 fprintf('%s(%3u)  ',buff{ii},les(ii));
end % ii
fprintf('\n');
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function indices = mystrfind( cellarr,le )
indices = [];
for ii = 1:length(cellarr),
 if strcmpi(cellarr{ii}(1),le),indices = [indices ii]; end
end % ii
end % eofunc
