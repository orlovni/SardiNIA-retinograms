%
% run_chop_up_sigMatrix_noClasses( 'CTmyRegistration_CLAHE_wholeTorso_Males_2gr',.5,'hehe' );
% run_chop_up_sigMatrix_noClasses( 'CTmyReg_rawscans_wTorso_females',.5,'hehe' );
% chop_up_AgeUniform_sigMatrix_noClasses( 'CTmyRegistration_CLAHE_STs_Males' );
%
function chop_up_AgeUniform_sigMatrix_noClasses( matname,TRportion,Ntiles,varargin ),

load( matname );
desired_classes = [];

info{2} = dataset_name;
info{3} = image_ids;
info{4} = image_paths;
info{5} = signature_labels;
info{6} = signature_matrix;
info{7} = desired_classes;

%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
% Out of given 'image_paths', get 1) 'whole_BLSA_list', 2) list with unique
% BLSA IDs ('uIDlist'), 3) 'nrep_uID': number of repetitions of 'uIDlist'
% in 'whole_BLSA_list'
% Information about unique BLSA IDs is a key in forming training set (that
% should be solely based on non-repetitive patient-visit entries)
[ uIDlist,nrep_uID,whole_BLSA_list,whole_Ages_list ] = get_unique_BLSAid_list( image_paths );

%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
% Basic stats on the data set
totalIDs = length( uIDlist );
fprintf('\nDataset: %s\n',matname);
fprintf('Total  number of BLSA IDs in the set:  %u\n',length(whole_BLSA_list));
fprintf('Number of unique BLSA IDs in the set:  %u\n',totalIDs);
fprintf('Portion intended for Training purpose: %.2f\n',TRportion);


info{8}  = uIDlist;
info{9}  = nrep_uID;
info{10} = whole_BLSA_list;
info{11} = whole_Ages_list;

%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
% I could save the TE (test) data as a single file (default option) using 
% the same path as for TR (train) data, OR save TE data as separate files.
% This would happen when there is a specified 4th parameter, the path to
% the TE data
if nargin > 3, 
 separate_exp_files = varargin{1}; 
 info{12} = separate_exp_files;
end


%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
% Refine desired training size
if abs(TRportion) < 1,desired_trSz = floor( abs(TRportion) * size(signature_matrix,2) ); end

if desired_trSz > size(signature_matrix,2),
 desired_trSz = floor( .95 * size(signature_matrix,2) );
end % desired_trSz

%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
Ncla = length(category_names);
existing_classes = 1:Ncla;
class_vec = signature_matrix(end,:);

%category_names = [];
info{1} = category_names;

%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
split_TrTest( TRportion,matname,info,Ntiles );

fprintf('ky-ky\n');
end % eofunc


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function [uIDlist,nrep,whole_BLSA_IDlist,whole_Ages_list] = get_unique_BLSAid_list( image_paths ),

% whole BLSA ID list:
n = length(image_paths); 
whole_BLSA_IDlist = [];
whole_Ages_list = [];
for ii = 1:n,
 [foldr,nm,ext] = fileparts( image_paths{ii} );
 [blsaID,age] = BLSAid_from_name( nm );
 
 whole_BLSA_IDlist = [ whole_BLSA_IDlist  blsaID ];
 whole_Ages_list   = [ whole_Ages_list    age ];
end % ii

% unique list:
uIDlist = unique( whole_BLSA_IDlist );
nu = length( uIDlist );

% number of repetitions:
nrep = zeros(1,nu);
for ii = 1:nu,
 nrep(ii) = length( find( uIDlist(ii) == whole_BLSA_IDlist ) );
end % ii
end % eofunc


% Extract numerical part of BLSA ID given the [tif] file name
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function [blsaID,age] = BLSAid_from_name( filename ),
pos1 = findstr( filename,'-' ); pos1 = pos1(1);
%pos2 = findstr( filename,'_' ); pos2 = pos2(1);
pos0 = findstr( filename,'_' ); pos2 = pos0(1); 

pos_end = findstr(filename,'_scanA');
if isempty(pos_end), pos3 = pos0(end);
else, pos3 = pos0(end-1);
end % new name convention (without _scanA at the end)

blsaID = str2num( filename(pos1+1:pos2-1) );

age = str2num( filename(pos3+1:pos_end-1) );
if ~isempty(age),return;end

age = str2num( filename(pos3+1:end) );
end % eofunc
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% This function was revised (Aug, 2013) to adopt changes in naming 
% covention for Registered images 
function [ID,visit] = get_ID_Visit( ci_name ),visit = []; ID = [];

% get ID:
IDstart = length('BLSA-') + 1;
IDend   = findstr( ci_name,'_' ); IDend = IDend(1) - 1;
ID = ci_name(IDstart:IDend);
ID = str2num(ID); 

% List of abbreviations meant to be as "Visit". This is because the files
% in this acquisition were renamed manually by Sokratis and some of the
% files were misspelled
visitAbbrevs = {'_Visit-','_Vist-'};

% get Visit:
for ii = 1:length(visitAbbrevs),
 VISstart = findstr( ci_name,visitAbbrevs{ii} );
 if isempty(VISstart),continue,end
 VISend = findstr( ci_name,'_' ); VISend = VISend(2) - 1;
 visit  = ci_name(VISstart:VISend); break;
end % ii
if isempty(visit),return; end % empty visit
visit = str2num(visit); 
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function desired = refine_desired_training_size( desired_trSz,Ntiles ),
desired = Ntiles * floor( desired_trSz/Ntiles );
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function [into_training,into_testing] = split_by_randomz( TRportion,whole_BLSA_list,whole_Ages_list )
% 1 thing: randomize whole_BLSA_list
neworder = randperm( length(whole_BLSA_list) );
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function [into_training,into_testing] = split_by_agesampling( TRportion,whole_BLSA_list,whole_Ages_list )
ages = unique( whole_Ages_list );
into_training = []; into_testing = [];

for agebins = ages,
 indexlist = find( agebins == whole_Ages_list );
 if isempty(indexlist),continue,end
%fprintf('Age: %u\n',agebins);
 
 found_ids = whole_BLSA_list( indexlist );
 n = length(found_ids);
 if isempty( found_ids ),continue,end
 if n == 1,
 %into_training = [into_training indexlist];
  candidatesINX_TR = indexlist(1);
  candidatesINX_TE = [];
 end % only one found
 if n >= 2,
  ntr = floor( n*TRportion ); if ntr < 1,ntr = 1;end
  nte = n - ntr;
  candidatesINX_TR = indexlist(1:ntr);
  candidatesINX_TE = indexlist(ntr+1:end);
 end % > 1
  
 candidatesIDs_TR = whole_BLSA_list( candidatesINX_TR );
 candidatesIDs_TE = whole_BLSA_list( candidatesINX_TE );
    
 % Overlaps between 'current' and 'cumulative' IDs.
  
 % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 % #1a: Should there be any overlap between newfound (local) TE_IDs -- 
 % unique(candidatesIDs_TE) -- and cumulative TR_IDs, unique(
 % whole_BLSA_list(into_training) ), then move those overlapped samples 
 % into local TR set
 komulative_TR_thing = whole_BLSA_list(into_training);
 overlapped_IDs = intersect( unique(candidatesIDs_TE),unique(komulative_TR_thing) );
 if ~isempty( overlapped_IDs ),
   
  % These indeces ... 
  indeces_4overlapped = get_indices4overlap( candidatesIDs_TE,overlapped_IDs );
  %%indeces_4overlapped = find(ismember( candidatesIDs_TE,overlapped_IDs ));
   
  indeces_4overlapped = candidatesINX_TE(indeces_4overlapped);
   
  % ...are to be removed from current test candidates...
   candidatesINX_TE = setdiff( candidatesINX_TE,indeces_4overlapped );
   
  % ...and added toward accumulated TR
  candidatesINX_TR = [candidatesINX_TR indeces_4overlapped];
   
  % Update IDs accordingly...
  candidatesIDs_TR = whole_BLSA_list( candidatesINX_TR );
  candidatesIDs_TE = whole_BLSA_list( candidatesINX_TE );
   
  clear komulative_TR_thing overlapped_IDs;
 end % overlapYES 1
  
  
 % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 % #1b: Should there be any overlap between local TR_IDs -- 
 % unique(candidatesIDs_TR) -- and cumulative TE_IDs, unique(
 % whole_BLSA_list(into_testing) ), then move those overlapped samples 
 % into local TE set.
 komulative_TE_thing = whole_BLSA_list(into_testing);
 overlapped_IDs = intersect( unique(candidatesIDs_TR),unique(komulative_TE_thing) );
 if ~isempty( overlapped_IDs ),
   
  indeces_4overlapped = get_indices4overlap( candidatesIDs_TR,overlapped_IDs );
  indeces_4overlapped = candidatesINX_TR(indeces_4overlapped);
   
  candidatesINX_TR = setdiff( candidatesINX_TR,indeces_4overlapped );
  candidatesINX_TE = [candidatesINX_TE indeces_4overlapped];
   
  clear komulative_TE_thing overlapped_IDs;
 end % overlapYES 1


 % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 % #2: Overlaps between current TE/TR IDs
 [overlapYES,candidatesINX_TR,candidatesINX_TE] = check_overlaps( candidatesIDs_TR,candidatesIDs_TE,candidatesINX_TR,candidatesINX_TE );
 if overlapYES,
  candidatesIDs_TR = whole_BLSA_list( candidatesINX_TR );
  candidatesIDs_TE = whole_BLSA_list( candidatesINX_TE );
 end % overlapYES 2
  
  % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 into_training = [into_training candidatesINX_TR];
 into_testing  = [into_testing  candidatesINX_TE];
  
end % agebins

candidatesIDs_TR = whole_BLSA_list(into_training);
candidatesIDs_TE = whole_BLSA_list(into_testing);
uTR = unique(candidatesIDs_TR);
uTE = unique(candidatesIDs_TE); 

end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function indeces_4overlapped = get_indices4overlap( candidatesIDs_TE,overlapped_IDs )
indeces_4overlapped = [];
for id = overlapped_IDs,
 if any( ismember( candidatesIDs_TE,id ) ),
  indeces_4overlapped = [indeces_4overlapped find(id==candidatesIDs_TE)];
 end % ismember
end % id
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function [overlapYES,candINX_TR,candINX_TE] = check_overlaps( candidatesIDs_TR,candidatesIDs_TE,candidatesINX_TR,candidatesINX_TE ),
candINX_TR = candidatesINX_TR;
candINX_TE = candidatesINX_TE;
uniqueID_TR = unique( candidatesIDs_TR );
uniqueID_TE = unique( candidatesIDs_TE );
overlap     = intersect( uniqueID_TR,uniqueID_TE );
overlapIND  = candidatesINX_TE (find(ismember( candidatesIDs_TE,overlap )));
overlapYES = ~isempty(overlap);
if overlapYES,
 list2move_TR = find(candidatesIDs_TE==overlap);
 list2move_TR = candidatesINX_TE(list2move_TR);
 candINX_TR = [candidatesINX_TR list2move_TR];
 candINX_TE = setdiff(candidatesINX_TE,list2move_TR);
end % overlap
end % eofunc


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function split_TrTest( TRportion,matname,info,Ntiles ),

matr = info{6};
desired_classes = info{7};

uc = unique( desired_classes );

new_class_vec = []; 
new_class_vecTr = []; 
new_class_vecTe = []; 
new_class = 0;

train_indx = [];
test_indx  = [];

%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
uIDlist  = info{8};
nrep_uID = info{9};
whole_BLSA_list = info{10};
whole_Ages_list = info{11};

[into_training,into_testing] = split_by_randomz( TRportion,whole_BLSA_list,whole_Ages_list );
%[into_training,into_testing] = split_by_agesampling( TRportion,whole_BLSA_list,whole_Ages_list );


%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
%IDs_to_train  = uIDlist( nrep_uID <= UPPERrepeatsLIM );
%IDs_to_exerim = uIDlist( nrep_uID >  UPPERrepeatsLIM );
% Get indices for those IDs:
%into_training = find_instances( whole_BLSA_list,IDs_to_train, 1:length(whole_BLSA_list) );
%into_testing  = find_instances( whole_BLSA_list,IDs_to_exerim,1:length(whole_BLSA_list) );

% ~ ~ debugging ~ ~ ~ ~
IDs_tr_currentClass = whole_BLSA_list( into_training );
IDs_te_currentClass = whole_BLSA_list( into_testing );

fprintf('\nTR set full size:\t%u',length(IDs_tr_currentClass));
fprintf('\nTE set full size:\t%u\n\n',length(IDs_te_currentClass));
 
print_1D_array( unique(IDs_tr_currentClass),15,sprintf('To training [unique IDs]. (Class: %u)',0) );
print_1D_array( unique(IDs_te_currentClass),15,sprintf('To testing  [unique IDs]. (Class: %u)',0) );
% ~ ~ debugging ~ ~ ~ ~


%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
category_names   = info{1};
dataset_name     = info{2};
imIDs            = info{3};
imPA             = info{4};
signature_labels = info{5};

train_indx = into_training;

%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
% form training data
form_train_data( matname,matr,new_class_vecTr,imIDs,imPA,train_indx,'train',signature_labels,category_names,dataset_name );
fprintf(': Train data saved :\n');

%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
% form Experiment/Prediction data
test_indx = into_testing;

% To have all TE data in one file, or in multiple per-image files?
if length( info ) >= 12,
 separate_exp_files = info{12};
 form_experiment_data( matname,matr,new_class_vecTe,imIDs,imPA,test_indx,signature_labels,category_names,dataset_name,separate_exp_files );
else,
 form_experiment_data( matname,matr,new_class_vecTe,imIDs,imPA,test_indx,signature_labels,category_names,dataset_name );
end % length >= 12

%form_experiment_data( matname,matr,new_class_vecTe,imIDs,imPA,test_indx,signature_labels,category_names,dataset_name );
fprintf(': Experiment data saved :\n');
end % eofunc


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function print_1D_array( arr,n,header ),
fprintf('\n%s\n',header);
narr = length(arr);
fprintf(':%u:\n',narr);
fprintf(' ##\t%s\n\n',num2str(1:n,'(%02u) '));
%fprintf('#\t%s\n',num2str(1:n,'%4u '));

cnt = 0; i0 = 1;
for ii = 1:n:narr,
 cnt = cnt + 1;
 i0 = ii; 
 i1 = i0 + n - 1;
 i1 = min(i1,narr);
 
%w/o numbering:
%fprintf('%s\n',num2str(arr(i0:i1),' %04u'));

 fprintf('(%02u)\t%s\n',cnt,num2str(arr(i0:i1),' %04u'));
end % ii
fprintf('\n');
end % eofunc


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function found = find_instances( multiple_instances,pattern,indices ), 
found = [];
for ii = 1:length( pattern ),
 incs = find( pattern(ii) == multiple_instances );
 found = [found indices(incs)];
end % ii
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function form_experiment_data( matname,matr,new_class_vec,imIDs,imPA,target_indx,signature_labels,category_names,dataset_name,varargin ),
signature_matrix = matr(:,target_indx);

image_ids   = imIDs( target_indx );
image_paths = imPA ( target_indx );

% sort out experiments according to age:
ages = signature_matrix(end,:);
[ages,indx] = sort( ages );
signature_matrix = signature_matrix(:,indx);
image_ids   = image_ids  (indx);
image_paths = image_paths(indx);


if nargin < 10,
 % DEBUG mode: save all experiments in one big chunk:
 savename = [matname '_test'];
 %%save(savename,'category_names','dataset_name','signature_labels','signature_matrix','image_ids','image_paths');
 save(savename,'category_names','dataset_name','signature_labels','signature_matrix','image_paths');

else,
    
 separate_exp_files = varargin{1};
 
 if separate_exp_files,
  experiment_folder = [matname '_train'];
  stream_experiments( image_paths,image_ids,signature_labels,signature_matrix,category_names,dataset_name,experiment_folder );
 else
  % roll back to 1st choice: a single TE file
  savename = [matname '_test'];
  save(savename,'category_names','dataset_name','signature_labels','signature_matrix','image_paths');
 end % if separate_exp_files
 
end % nargin < 10 ?

% Go through the list of tile names (image_paths) and save experiments by
% the image name
%stream_experiments( image_paths,image_ids,signature_labels,signature_matrix,category_names,dataset_name );
%return;

end % eofunc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function mymkdir( outpath )
if ~exist( outpath ),
 mkdir( outpath );
end % ~exist
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function stream_experiments( image_paths0,image_ids0,signature_labels,signature_matrix0,category_names,dataset_name,experiment_folder ),
image_paths0 = image_paths0'; imagenames = {};
exp_cnt = 0;
current_batch = 0;
signature_matrix = [];
[nrows,ncols] = size( signature_matrix0 );
mymkdir( experiment_folder );

for ii = 1:length( image_paths0 ),
 [foo,curr] = fileparts( image_paths0{ii} );
 foo2 = strfind( imagenames,curr(1:end-4) );
 
 current_batch = current_batch + 1;
 signature_matrix = signature_matrix0( :,ii );
 image_ids = image_ids0(ii);
 image_paths = image_paths0{ii};
 exp_filename = fullfile( experiment_folder,curr );
 

%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
% setting up a new batch
 if length(foo2) == 0 | (length(foo2) >= 1 & length(foo2{end})==0),
     
%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
  % save pre-existing batch -- should such one exist...
  if ~isempty( signature_matrix ) %& (ncols >= Ntiles & rem(ncols,Ntiles)==0),
   save( exp_filename,'category_names','dataset_name','signature_labels','signature_matrix','image_ids','image_paths' );
  end
  continue;
  
  pti = strfind(curr,'_'); pti = pti(end)-1; 
  patt = curr(1:pti);
  
  imagenames = [imagenames; patt];
  exp_cnt = exp_cnt + 1;
  new_experiment = ['exp' num2str(exp_cnt,'%03u') '_'];
  
%~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~
% clean up space for the new batch
  current_batch = 1;
  image_paths = {};
  image_ids = [];
  signature_matrix = [];
  exp_filename = [new_experiment patt];
 
 end % if
 
end % ii
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function my_mkdir( foo ),
 if ~exist( foo,'dir' ),
  [status,message,messageid] = mkdir( foo );
  if ~logical(status),error(['** cannot make folder "' foo '" (' message ') **']);end
 end % exist
end % eofunc


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function imagenames = extract_imagenames_from_tilenames( image_paths,Ntiles ),
image_paths = image_paths'; imagenames = {};
for ii = 1:length( image_paths ),
 [foo,curr] = fileparts( image_paths{ii} );
 foo2 = strfind( imagenames,curr(1:end-4) );
 if length(foo2) == 0 | (length(foo2) >= 1 & length(foo2{end})==0),
  patt = curr(1:end-4);
  imagenames = [imagenames; patt];
 end % if
end % ii
end % eofunc

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function [indlist,image_paths2] = findalltiles( image_paths,patt ),
indlist = []; image_paths2 = {}; cnt = 0;
for ii = 1:length( image_paths ),
 k_tiles = findstr( image_paths{ii},patt );
 if ~isempty(k_tiles), indlist = [indlist ii]; 
 else, cnt = cnt+1; image_paths2{cnt} = image_paths{ii};
 end % ~isempty
end % ii
end % eofunc


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function mn = min_class_Size( perClass,uc ),
class_sizes = [];
for cla = uc,
 current_class = perClass(cla).index;
 n = length(current_class);
 class_sizes = [class_sizes n];
end % cla
mn = min(class_sizes);
end % eofunc


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
function form_train_data( matname,matr,new_class_vec,imIDs,imPA,target_indx,hint,signature_labels,category_names,dataset_name ),
savename = [matname '_' hint];

signature_matrix = matr(:,target_indx);

image_ids = imIDs(target_indx);
image_paths = imPA(target_indx);

ages = signature_matrix(end,:);
[ages,indx] = sort( ages );

signature_matrix = signature_matrix(:,indx);
image_ids   = image_ids  (indx);
image_paths = image_paths(indx);

%%save( savename,'category_names','dataset_name','signature_labels','signature_matrix','image_ids','image_paths' );
save( savename,'category_names','dataset_name','signature_labels','signature_matrix','image_paths' );
end % eofunc

