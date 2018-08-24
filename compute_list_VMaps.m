% Given name of folder with images and name of destination folder for
% signatures (where some signatures might be partially computed) generates
% corresponding file names for each list (source and destination)
%
% Example 1: ( using within matlab environment )
%
% curDir = '~/Documents/MATLAB/Fundus_images/_Comp_chain/';
% compute_list_VMaps([curDir filesep '../_extracted_skeletonized_CLAHE_0_999.log'],[curDir filesep '../_computed_sigs_0_999']);
% compute_list_VMaps([curDir filesep '../_extracted_skeletonized_29img'],[curDir filesep '../_computed_sigs_29img']);
%
% Example 2: ( using w/o matlab environment )
% matlab < run_compute_list_VMaps.m > aaa_Test_segmentation.log & 
%
% Edits: Oct 2016; Jun 2017

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function compute_list_VMaps (varargin)

    % Folder that contains image files to be processed
TaskFolder = varargin{1}; 
    % Folder that contains files that already processed
ReadyFolder = varargin{2};


dataEXT = '.png';


fprintf('\n\n')
fprintf('image files from: %s\n',TaskFolder);
fprintf('*.mat files to  : %s\n',ReadyFolder);
fprintf('\n\n')

if ~exist(ReadyFolder),
 smd = mkdir(ReadyFolder);
 fprintf('creating folder: status %u\n',smd);
else,disp([ReadyFolder ' already exist:: using existing folder']);
end

    % Containing of both folders (includes dirnames, other files, etc.)
task_files_all  = dir(fullfile(TaskFolder,['*' dataEXT])); 
ready_files_all = dir(fullfile(ReadyFolder,'*.mat')); 

fprintf('estimate size of both folders:\n');
whos task_files_all ready_files_all


    % Cleaning file lists (exclude dirnames, other files, etc.)
%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Current segm/skel maps:
[task_files,sext]  = cleaning_flist(task_files_all,dataEXT);

%[task_files,sext]  = cleaning_flist(task_files_all,'.png');

% tiffs:
%[task_files,sext]  = cleaning_flist(task_files_all,'.tiff');

%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

[ready_files,dext] = cleaning_flist(ready_files_all,'.mat');

%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
if length(task_files)<1, error(':: no source files, abort mission ::');end
if (length(task_files)>=1) & (isempty(task_files{1})), error(':: no source files, abort mission ::');end

fprintf('estimate size after cleaning:\n');
whos task_files ready_files
%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

if prod(size(ready_files))==1 & isempty(ready_files{1}),
 [files2process,files2save] = generate_diff_list(task_files,{},TaskFolder,ReadyFolder,sext,dext);
else,
 [files2process,files2save] = generate_diff_list(task_files,ready_files,TaskFolder,ReadyFolder,sext,dext);
end % if
fprintf('.. list has been generated ..\n\n');
fprintf('   files to process: %u\n',length(files2process));
if length(files2process) < 1, error('no files to process in src folder: abort'); end

% warning off MATLAB:divideByZero
DDB = getSigsA( files2process,ReadyFolder,dataEXT );

end % eofun


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [cleaned_list,ext] = cleaning_flist(flist,fext1),
ext = [];
NN = length(flist); cleaned_list = cell(1);
for ii = 1:NN,
 FName = flist(ii).name; dirFlag  = flist(ii).isdir;
 [dr,fn,ext] = fileparts(FName);
 
%fn = lower(fn); ext = lower(ext);
 if (strcmpi(FName,'.DS_Store') | ~strcmpi(ext,fext1) | dirFlag), continue; end
 if isempty(cleaned_list{1}),
  cleaned_list{1} = fn;
 else,
  cleaned_list = [cleaned_list cellstr(fn)];
 end
end % for
%ext = fext1;
end % eofun

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [files2process,files2save] = generate_diff_list(task_files,ready_files,SourceFold,DestFold,sext,dext),
N_task  = length(task_files);
N_ready = length(ready_files);
files2process = {[]}; files2save = {[]};
for itask = 1:N_task,
% Compare current task_file against list of computed (ready_files)   

 if isempty(strmatch(lower(task_files{itask}), lower(ready_files))),
  if isempty(files2process{1}), 
   files2process{1} = fullfile(SourceFold,[task_files{itask} sext]);
   files2save{1}    = [task_files{itask} dext];
  else, 
   files2process = [files2process cellstr(fullfile(SourceFold,[task_files{itask} sext]))];
   files2save    = [files2process cellstr([task_files{itask} dext])];
  end % if isempty2
 end % if isempty
end % for
end % eofun



% compile name for no-vessel dump:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function no_vessels_name = get_no_vessels_name( foo ),
[pa,na] = fileparts( foo );
k_filesep = strfind( pa,'/' );
k_filesep = 1 + k_filesep(end);
foo = foo( k_filesep:end );
k_dot = strfind( foo,'.' );
if ~isempty(k_dot),foo = foo(1:(k_dot-1));end
%no_vessels_name = ['no_vessels' foo '.mat'];

% N: June 19,2017
[ff,fn] = fileparts( foo );
%%no_vessels_name = ['no_vessels' foo '.mat'];
no_vessels_name = ['no_vessels' ff '.mat'];
end % eofun


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function ID = get_imageID( nm ),
[foo,na,ext] = fileparts( nm );

% There are two types of name convention: 
% IM0001.tif (~testType) and 0001_test.tif (testType)

k0 = strfind( na,'IM' );
k1 = strfind( na,'_' );
s = na((k0+length('IM')):(k1-1));
 
ID = str2num(s);
end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [DDB,siglabs] = getSigsA( imgNames,dataNameStamp,dataEXT )
%   INPUT:          imgNames    - filenames of input images.  i.e. 'image1.tif'
%   OUTPUT:         DDB          - matrix of signatures.  columns are images, rows are signatures.
%                                   so, DDB(1,1) is the first image's first signature and DDB(5,1) is
%                                   the first image's fifth signature.  DDB(1,5) is the fifth image's
%                                   first signature
%                   siglabs     - names of all of the signatures outputed in DDB
%   Written by Nikita Orlov - 2003-2005

[imgSize imgSize] = size(imread(imgNames{1}));
img_quantity      = length(imgNames);         % how many images

if ~exist(dataNameStamp),
 stmkdir = mkdir(dataNameStamp);
 
if ~stmkdir,
 if7 = exist(fullfile(pwd,dataNameStamp));
 matfiles_computed = {};
 if if7 ~= 7, fprintf('fld :: %s\n',dataNameStamp); error('!! can''t create folder; folder does''t exist !!'); end
end % if ~stmkdir

else,
 % what is in there? what IDs in there?
 
 dtaFolder = fileparts(imgNames{1});
 matfiles_computed = dir( fullfile( dtaFolder,'*.mat' ) );
%matfiles_computed = dir( fullfile( dataNameStamp,'*.mat' ) );
 matfiles_computed = {matfiles_computed.name};
 nm = length(matfiles_computed);
 matfiles_computedIDs = [];
 for ii = 1:nm,
  ID = get_imageID( matfiles_computed{ii} );
  matfiles_computedIDs = [matfiles_computedIDs ID];
 end % ii
 
end % if ~exist(dataNameStamp)

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
foo = imgNames{1};

% compile name for no-vessel dump:
no_vessels_name = get_no_vessels_name( foo );

foo = fileparts(foo);

% Compute number of segm masks & mat files in the data folder
[n_segm,n_mat] = get_numberOf_segmNmats( foo );

% Do I have all mean calibers deposited to the data folder?
% if n_mat == n_segm I have matching number of files and thus I can read in
% all mean calibers and compute their mean as Uber-mean (eg masterCaliber)
% for the whole set. This is what I needed.
%
% Otherwise, I proceed with the 'preliminary analysis' step, followed again
% by reading in all mean calibers and computing their mean as Uber-mean.
all_mean_calibers_in_place = (n_mat == n_segm);

% The above criterion is good only when all skeleton images were of
% acceptable quality. Iff some skeletons were skipped (on the grounds of
% poor quality), then n_mat < n_segm. This is a real life scenario, and
% that's why I use a FLAG 'calibers_completed' as a new criterion.
%%%load no_vessels;
if exist( no_vessels_name,'file' ),
 try,
  load( no_vessels_name );
 catch, 
  fprintf('>> file %s did not load <<\n\n',no_vessels_name);
  calibers_completed = true;
 end % try
end

if exist( 'calibers_completed','var' ),
 all_mean_calibers_in_place = calibers_completed;
 load master_caliber;
else, 
 all_mean_calibers_in_place = false;
 noVessels_list_skel = {};
 noVessels_list_segm = {};
 yesVessels_list_segm = {};
 calibers_completed = false;
 cnt = 0; 
 noVessels_cnt  = 0;
 yesVessels_cnt = 0;
end % if exist


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% This preliminary analysis takes all images to analyse their calibers
if all_mean_calibers_in_place & ~exist( 'masterCaliber','var' ),
    
 masterCaliber = get_masterCaliber( foo );
 
elseif ~exist( 'masterCaliber','var' ),
    
 masterCaliber = [];
 
% these values have been read from the existing file ( 'no_vessels_name' )
 
 for ii = 1:img_quantity, cnt = cnt + 1;  
     
    if ii == img_quantity, calibers_completed = true; end
         
   %if ii < 1269,continue;end
   %if ii < 6,continue;end
   %if cnt ~= 35,continue;end
    
    tic0 = cputime;
    currentImageName = imgNames{ii};
    [isSkelName,skeletonImgName] = analyseIMname( currentImageName );
    if isSkelName,continue;end
    
    
    if exist( no_vessels_name,'file' ),
        
    % Check: if the corresponding mat file already in place and thus I am
    % to skip the image file
     if exist( 'yesVessels_list_segm','var' ),
       skip_existing_mat = skip_mat_YN( matfiles_computedIDs,currentImageName );
     else, skip_existing_mat = false;
     end % variable yesVessels_list_segm exists
     
    else, skip_existing_mat = false;    
    end % no_vessel_name exists

    
    segmentedImgName = currentImageName;
    fprintf('image %03u / %03u [%s]',ii,img_quantity,currentImageName);
    isaborted = false;
    
    if ~skip_existing_mat,
    % Does not compute sigs. Only saves average map branch thickness to the data folder 
     [DDBbufer,sigSizes,isaborted] = sigchain_VMaps_segm( segmentedImgName,skeletonImgName,masterCaliber );
     tmStamp = getEtime(cputime-tic0);
    else,
     tmStamp = 'mat exists,image skipped';
    end %
    
    % show the file name for AVG caliber just saved ( or aborted )
    meanCaliberFname = gen_AVGthick_fname( segmentedImgName,isaborted );
    
    if isaborted,
      tmStamp = 'noVessels,image skipped';
      noVessels_cnt = noVessels_cnt + 1;
      noVessels_list_skel{noVessels_cnt} = skeletonImgName;
      noVessels_list_segm{noVessels_cnt} = segmentedImgName;
 %    save('no_vessels','noVessels_list_skel','noVessels_list_segm','calibers_completed');
    else,
      yesVessels_cnt = yesVessels_cnt + 1;
      yesVessels_list_segm{yesVessels_cnt} = segmentedImgName;
    end % if isaborted
    
    save(no_vessels_name,'noVessels_list_skel','noVessels_list_segm','yesVessels_list_segm','calibers_completed','noVessels_cnt','yesVessels_cnt','cnt');
    
    fprintf('  (%s )\n',tmStamp);
    
 end % ii
 
 masterCaliber = get_masterCaliber( foo );
end % if all_mean_calibers_in_place

%fprintf( '\n\nMasterCaliber over %u images has been computed\n\n',img_quantity );
fprintf( '\n\nMasterCaliber over %u images has been computed\n\n',n_mat );

TSSstructure = get_thick_segm_skel( foo,dataEXT );

[ff,nn] = fileparts( foo );

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% This ad hoc chain operates when image calibers have been already written
cnt = 0;
very_1st_sample_computed = true;

nTSS = length( TSSstructure );
img_quantity = nTSS;

for ii = 1:nTSS, %%img_quantity,   
    tic0 = cputime;
    
    segmentedImgName = TSSstructure(ii).segm;
    skeletonImgName  = TSSstructure(ii).skel;
    thicknessMATname = TSSstructure(ii).thick;
    
    % reading 'avgThick':
    try,
%    curr = load( fullfile( foo,thicknessMATname ));
     load( fullfile( foo,thicknessMATname ));
    catch
     fprintf('foo: %s\n',foo);
%    fprintf('ff: %s\n',ff);
%    fprintf('nn: %s\n',nn);
    end % try
    
    if 1 ==0
    currentImageName = imgNames{ii};
    % Analyse image name: segmented or skeleton? 
    [isSkelName,skeletonImgName] = analyseIMname( currentImageName );
    % if currentImageName is a skeleton image, skip 
    if isSkelName,continue;end
    segmentedImgName = currentImageName;
    end % 1==0
    
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
    fprintf('image %03u / %03u [%s]',ii,img_quantity,segmentedImgName);
   %if ismember( skeletonImgName,noVessels_list_skel ),
    if avgThick < 1e-2,
      skipImage = logical(1); fprintf('  ( skipped )\n'); continue;
    end %
    
    [DDBbufer,sigSizes] = sigchain_VMaps_segm( fullfile(foo,segmentedImgName),fullfile(foo,skeletonImgName),masterCaliber );
    
    % NOT SURE ABOUT THIS BLOCK ...
    if very_1st_sample_computed, % was: if ii == 1
     [mddb,nddb] = size(DDBbufer); % memorize the right size of signature vec.
     siglabs = get_siglabs( sigSizes ); 
     labels_memo = gen_long_labels;
     very_1st_sample_computed = false;
    else,
     % if problem with sig vector:   
     if size(DDBbufer,1) ~= mddb, fprintf('  ( skipped )\n'); continue;end
    end

    cnt = cnt + 1;   
    DDB(:,cnt) = DDBbufer(:);
    tmStamp = getEtime(cputime-tic0);
    fprintf('  (%s )\n',tmStamp);
    
% dump computed file as *.mat
%   dummyNm = imgNames{ii};
    dummyNm = segmentedImgName;
    [bufDir,bufFname,bufExt] = fileparts(dummyNm); 
    dummyDDB = DDB(:,cnt);
    MATname = fullfile(dataNameStamp,[bufFname '.mat']);
    save(MATname,'dummyDDB','siglabs','dummyNm','labels_memo');
end % ii
end % eofun



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function TSSstructure = get_thick_segm_skel( fldName,dataEXT ),
thickNameLST = get_name_lst( fullfile( fldName,'*AVGthick*.mat' ) );
segmNameLST  = get_name_lst( fullfile( fldName,['*segm*' dataEXT] ) );
skelNameLST  = get_name_lst( fullfile( fldName,['*skel*' dataEXT] ) );
TSSstructure = struct( 'thick',[],'segm',[],'skel',[] );

% The shortest list is always '*.mat':
n = length( thickNameLST );
cnt = 0;
for ii = 1:n,
 currentMAT = thickNameLST{ii};
 
 % Here is the typical mat file structure:
 %  51011-0-Age_76-Female_AVGthick__lmse_scales2.mat
 % Get the prefix, the identifying data
 IDstamp = get_IDprefix( currentMAT );
 
 % Detect matching segm and skel files for the IDstamp:
 matchingSEGM = get_img_match( segmNameLST,IDstamp );
 matchingSKEL = get_img_match( skelNameLST,IDstamp );
 
 if ~isempty( matchingSEGM ) & ~isempty( matchingSKEL ),
  cnt = cnt + 1;
  TSSstructure(cnt).thick = currentMAT;
  TSSstructure(cnt).segm  = matchingSEGM;
  TSSstructure(cnt).skel  = matchingSKEL;
 end % ~isempty
end % ii

end % eofun

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function nameLST = get_name_lst( pattn ),
foo = dir( pattn ); nameLST = { foo.name };
end % eofun

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function IDstamp = get_IDprefix( fName )
k = strfind( fName,'_' ); k = k(2)-1;
IDstamp = fName( 1:k );
end % eofun

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function matchFILE = get_img_match( nameLST,IDstamp )
k = strfind( nameLST,IDstamp );
k = find( ~cellfun(@isempty,k) );
if isempty(k),matchFILE = [];return;end
if length(k) > 1,se=['>> get_img_match: ' IDstamp ': multiple matches,abort <<'];error(se);end
matchFILE = nameLST{k};
end % eofun


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function skip_existing_mat = skip_mat_YN( matfiles_computed_IDs,currentImageName ),

skip_existing_mat = false;

[pa,na,ex] = fileparts( currentImageName );

% current ID:
ID = get_imageID( na );

if ismember( ID,matfiles_computed_IDs ),
 skip_existing_mat = true;
end

end % eofun


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function meanCaliberFname = gen_AVGthick_fname( mask_path,isaborted ),
SG = '_segm_';
[fileDir,fileNm,fileExt] = fileparts( mask_path );
k0 = strfind( fileNm,SG ) - 1;
k1 = k0 + length(SG);
part1 = fileNm(1:k0);
part2 = fileNm(k1:end);

if isaborted,
 meanCaliberFname = [part1 '_AVGthick_' part2 '_SKIPPED'];
else,
 meanCaliberFname = [part1 '_AVGthick_' part2];
end % if isaborted

meanCaliberFname = fullfile( fileDir,meanCaliberFname );
end % end function


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function masterCaliber = get_masterCaliber( path_indivMats ),
foo = dir( fullfile(path_indivMats,'*.mat') ); n = length(foo);
thick_array = [];
for ii = 1:n,
 fn = foo(ii).name;
 load(fullfile(path_indivMats,fn));
 thick_array = [thick_array avgThick];
end % ii
masterCaliber = mean( thick_array );
save master_caliber masterCaliber thick_array
end % eofun


% Compute number of segm masks & mat files in the data folder
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [n_segm,n_mat] = get_numberOf_segmNmats( foo ),
%strcomm = ['ls ' foo '/*.mat | wc']; [sta_mat,n_mat] = system( strcomm );

[ff,nn] = fileparts(foo); strcomm = ['ls ' nn '/*.mat | wc']; [sta_mat,n_mat] = system( strcomm );

if strcmpi(n_mat(1:3),'ls:') | strcmpi(n_mat(2:4),'ls:'),
 n_mat = 0;
else,
 if ischar( n_mat ),n_mat = str2num( n_mat );end
 n_mat = n_mat(1);
end

%if ~isempty(n_mat), n_mat = str2num( n_mat ); n_mat = n_mat(1);
%else, n_mat = 0;
%end

%strcomm = ['ls ' foo '/*segm* | wc'];
strcomm = ['ls ' foo '/*skel*.png | wc'];
strcomm = ['ls ' nn '/*skel*.png | wc'];
[sta_segm,n_segm] = system( strcomm );

if strcmpi(n_segm(1:3),'ls:') | strcmpi(n_segm(2:4),'ls:'),
 n_segm = 0;
else,
 n_segm = str2num( n_segm );
 n_segm = n_segm(1);
end

%if ~isempty(n_segm), n_segm = str2num( n_segm );n_segm = n_segm(1);
%else, n_segm = 0;
%end
end % eofun


% Analyse image name: segmented or skeleton?
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [isSkelName,skeletonImgName] = analyseIMname( currentImageName )
SK = '_skel_';
ske = strfind( currentImageName,SK );
if isempty( ske ),
 skeletonImgName = [];
 isSkelName = logical(0);
 % compile skeleton image name:
 SG = '_segm_';
 k = strfind( currentImageName,SG );
 part1 = currentImageName( 1:(k-1) );
 part2 = currentImageName( (k+length(SG)):end );
 skeletonImgName = [part1 SK part2];
else,
 skeletonImgName = currentImageName;
 isSkelName = logical(1);
end % isempty
end % eofun

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function im2 = randomlyChangeSize(img),
image_scales = [0.5:0.025:1.3 1.5:.2:2];
rand('state', sum(100*clock));
image_scales = image_scales(randperm(length(image_scales)));
image_scales = image_scales(1);
im2 = imresize(img,image_scales,'bicubic');
end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [st,hh,mm,ss]=getEtime(tm),
hh = floor(tm/3600); hrem = rem(tm,3600);
mm = floor(hrem/60); ss = rem(hrem,60); 
st = sprintf('%02i:%02i:%02i.%02i',hh,mm,fix(ss),round(1e2*(ss-fix(ss))));
end % eof


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Creating all of the labels for all of the signatures.  
function siglabs = get_siglabs( sigSizes ),
%function siglabs = get_siglabs(sigSizes,ch_),
uniqueLabNames = {...
	['wholeM_' 'Tortuosity15 '], ...               
	['wholeM_' 'FractalDim08 '], ...               
	['wholeM_' 'JunctionAnalysis03 '], ...               
	['thin_' 'Tortuosity15 '], ...               
	['thin_' 'FractalDim08 '], ...               
	['thin_' 'JunctionAnalysis03 '], ...               
	['thick_' 'Tortuosity15 '], ...               
	['thick_' 'FractalDim08 '], ...               
	['thick_' 'JunctionAnalysis03 '], ...               
};

siglabs = [];
for cnt = 1:length(sigSizes),
 siglabs = addSigLabel(siglabs,uniqueLabNames{cnt},sigSizes(cnt));
end
siglabs = siglabs';
end % eofun

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generates notes/comments for the feature labels  
function labels_memo = gen_long_labels,
labels_memo = {...
    'The three feature families are:',...
    ' Junction Analysis Features (vector of length 3)',...
    ' Tortuosity Features (TFs), (vector of length 15)',...
    ' Fractal Dimension Features (vector of length 8)',...
    '',...
    'Each feature family is computed on the whole vessel map,',... 
    ' and for every "vessel thickness map"',...
    '',...
    'Junction Analysis Features (JAFs):',...
    ' #1: terminal points ( T), # of',...
    ' #2: bifurcation points ( Y ), # of',...
    ' #3: crossing points ( X ), # of',...
    '',...
    'Tortuosity Features (TFs):',...
    ' #01: ratio of arcLength to cordLength, sampled midline',...
    ' #02: ratio of arcLength to cordLength, interpolated midline',...
    ' #03: total curvature, sampled midline',...
    ' #04: total curvature, interpolated midline',...
    ' #05: total squared curvature, sampled midline',...
    ' #06: total squared curvature, interpolated midline',...
    ' #07: total curvature / arc length, sampled midline',...
    ' #08: total curvature / arc length, interpolated midline',...
    ' #09: total squared curvature / arc length, sampled midline',...
    ' #10: total squared curvature / arc length, interpolated midline',...
    ' #11: total curvature / chord length, sampled midline',...
    ' #12: total curvature / chord length, interpolated midline',...
    ' #13: total squared curvature / chord length, sampled midline',...
    ' #14: total squared curvature / chord length, interpolated midline',...
    ' #15: ratio of sampled arc length to interpolated arc length',...
    '',...
    'Fractal Dimension Features (FDFs):',...
    ' Scale-specific box count measures are defined in terms of x[i],y[i],',...
    ' where x[i] := Resolutions[i], and y[i] = log(boxCounts[i])./log(x[i])',...
    ' ##1:8: y[i]',...
    };
return;
labels_memo{1} = 'Notes on feature labels';
labels_memo{2} = 'Junction Analysis Features (JAFs), ##1:3';
labels_memo{3} = ' #1: terminal points, # of';
labels_memo{4} = ' #2: bifurcation points ( Y ), # of';
labels_memo{5} = ' #3: crossing points ( X ), # of';
labels_memo{6} = ' JAFs: computed on whole vessel map & for every vessel thickness';
labels_memo{7} = 'Tortuosity Features (TFs), ##1:15';
labels_memo{8} = ' #01: ratio of arcLength to cordLength, sampled midline';
labels_memo{9} = ' #02: ratio of arcLength to cordLength, interpolated midline';
labels_memo{10} = ' #03: total curvature, sampled midline';
labels_memo{11} = ' #04: total curvature, interpolated midline';
labels_memo{12} = ' #05: total squared curvature, sampled midline';
labels_memo{13} = ' #06: total squared curvature, interpolated midline';
labels_memo{14} = ' #07: total curvature / arc length, sampled midline';
labels_memo{15} = ' #08: total curvature / arc length, interpolated midline';
labels_memo{16} = ' #09: total squared curvature / arc length, sampled midline';
labels_memo{17} = ' #10: total squared curvature / arc length, interpolated midline';
labels_memo{18} = ' #11: total curvature / chord length, sampled midline';
labels_memo{19} = ' #12: total curvature / chord length, interpolated midline';
labels_memo{20} = ' #13: total squared curvature / chord length, sampled midline';
labels_memo{21} = ' #14: total squared curvature / chord length, interpolated midline';
labels_memo{22} = ' #15: ratio of sampled arc length to interpolated arc length';
labels_memo{23} = ' TFs: computed on whole map  & for every vessel thickness';
labels_memo{24} = 'Fractal Dimension Features (FDFs), ##1:8';
labels_memo{25} = ' Resolutions = 1/BoxSize (that is 1,2,3,4,...); x[i] := Resolutions[i]';
labels_memo{26} = ' ##1:8: y[i] = log(boxCounts[i])./log(x[i])';
labels_memo{27} = ' FDFs: computed on whole map  & for every vessel thickness';
end % eofun

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function siglabs = addSigLabel(siglabs0,pattern,sigSize), siglabs = siglabs0;
buff = []; for ii = 1:sigSize, buff = [buff; [pattern sprintf('%2i',ii)]]; end
if ~isempty(buff), buff = cellstr(buff);
 else, buff = {pattern};
end
if isempty(siglabs), siglabs = buff;
else, siglabs(size(siglabs,1)+1:size(siglabs,1)+size(buff,1)) = buff;
end
end % eofun


function print_SigLab(n1,n2,sigLab),
for ii = n1:n2, fprintf('%s; ',sigLab{ii}); end
fprintf('\n\n');
end % eofun

