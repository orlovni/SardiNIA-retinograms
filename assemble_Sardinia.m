% ATTENTION! 
% 'assemble_Sardinia.m' uses different default rules on data
%
% New rules for MAT and TIFF folder structures:
% they contain only subfolders corresponding to class names, as following:
%   subject_01, subject_02, ...


%
% Examples: 
% assemble_Sardinia('../_computed_sigs_Sardinia6500','Fundus_Good1andBad2_groups');


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function assemble_Sardinia(varargin)
datapath = varargin{1};
if nargin > 1, output_mat_file_name = varargin{2}; else, output_mat_file_name = 'cla_du_DDB_file'; warning(['no output name -- use ' output_mat_file_name 'as default']); end
if nargin > 2, image_source_folder = varargin{3}; 
else, 
 image_source_folder = datapath;
 %error('no image_source_folder: abort'); 
end

output_mat_file_name = ['TR_' output_mat_file_name];

% Count all classes
disp(datapath);
[matFold_Names, Ncla]   = get_class_list(datapath);
if length(matFold_Names)==0,error(['Problem: ' datapath ' did not find data files']);end

disp(image_source_folder);
%[imgFold_Names, Ncla]   = get_class_list(image_source_folder);


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% This loop goes over the image source folder (to be containing subfolders with original image files).
% The goal here is to construct cell array ('orig_images'); its cells would be lists of images of the class.
im_IDs = uint16([]); cnt_IDs = 0; %cnt_IDs = 1000;
for cla = 1:Ncla,
 filebuf = dir(fullfile( datapath,matFold_Names{cla} ));
 
% Read all images from current folder... 
 last_image_index_string = [];
 %last_class_templ = [];
 current_class_index = 0;
 cnt_matfile = 0;
 for ii = 1:length(filebuf),
  fname = filebuf(ii).name; dirFlag = filebuf(ii).isdir; % image tiff file
  if dirFlag | (fname(1)=='.'), continue; end
  [foo1,nm,fext] = fileparts(fname);

  cnt_matfile = cnt_matfile + 1;
  current_image_index_string = num2str(cnt_matfile);

% Here we figure out what is the current_image_index and whether it has
% just changed ( since the last loop iteration ) 
  if isempty(last_image_index_string), 
      % The very 1st image in current class_folder: increase image_index ...
   last_image_index_string = current_image_index_string;
   current_class_index = cla;
   cnt_IDs = current_class_index; % cnt_IDs is the class index ...
   
   im_IDs = [im_IDs cnt_IDs];
   
      % Current_image_index ~= Last_image_index : the same effect -- increase image_index...
  elseif ( str2num(current_image_index_string) ~= str2num(last_image_index_string) )   
    
   last_image_index_string = current_image_index_string;
   im_IDs = [im_IDs cnt_IDs];
   
      % Current_image_index == Last_image_index : do nothing ...
  elseif ( str2num(current_image_index_string) == str2num(last_image_index_string) )     

   im_IDs = [im_IDs cnt_IDs];
  end %% if isempty...
   
  orig_images{cla}.name{cnt_matfile} = fname;
  orig_images{cla}.image_ids(cnt_matfile) = cnt_IDs;
 end % ii
end % cla



% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
DDB = []; trainImgs_train = {}; cnt = 0;
for cla = 1:Ncla, cnt_class = 0;

% This loop goes over the real data; each file would be compared with the structure 'orig_images'
% and will be assigned a name accordingly.
allfiles = dir(fullfile( datapath,matFold_Names{cla} )); N = length(allfiles);
 
class_vec = [];
for ii = 1:N,
  dummy = allfiles(ii).name; dirFlag  = allfiles(ii).isdir;
  [dr,fn,ext] = fileparts(dummy);
  if (strcmpi(dummy,'.DS_Store') | ~strcmpi(ext,'.mat') | dirFlag), continue; end
  load( fullfile(datapath, matFold_Names{cla}, dummy) );
  
  if isempty( dummyDDB ),continue;end
  
  siglabs = check_labels(siglabs);

  currClass = cla; 
  current_imID = [ orig_images{cla}.image_ids ];
  
  dummyDDB = convert_to_one_column(dummyDDB);
  cnt = cnt + 1;
  cnt_class = cnt_class + 1;
  im_path{cnt} = dummyNm; sig_labels = siglabs;
  image_ids(cnt) = current_imID(cnt_class); 
  try
  DDB = [DDB [dummyDDB; currClass]];
  catch
   fprintf('\n\n');
   whos DDB dummyDDB
  end % try
  fprintf('class: %02u %04u / %04u  (#%04u)\t%s\t\t>> {imID:%u} \n',currClass,ii,N,cnt,dummy,current_imID(cnt_class));
end % ii
fprintf('Class # %u: %u files\n',cla,cnt);

end % cla
class_vec = DDB(end,:);
image_ids = uint16(image_ids);

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        % Now I want to toss DDB columns so that classes would be grouped together
[DDB,class_vec,im_path] = tossDDB(DDB,class_vec,im_path);
 
trainImgs = im_path;
%save(output_mat_file_name,'DDB','sig_labels','trainImgs','class_vec','image_ids'); 

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% New saving style (the one implemented with 'transfer to correct vars' script)
category_names   = matFold_Names;
dataset_name     = output_mat_file_name;
signature_labels = sig_labels;
signature_matrix = DDB;
image_paths      = trainImgs;
notes = char(labels_memo(:));

%save(output_mat_file_name,'category_names','dataset_name','image_paths','signature_labels','signature_matrix','image_ids','notes');

% NB:
% Stop saving 'image_ids'!
% Providing 'image_ids' is dangerous, as it tends to faulter implementation
% of splits ( function 'updateSplitDirs.m' ). The function
% 'updateSplitDirs.m' expects variables 'image_ids' to be all different and
% to mimic different multi-tile images, while in reality the 'image_ids'
% presented had instead all the same values within the class (!).
save(output_mat_file_name,'category_names','dataset_name','image_paths','signature_labels','signature_matrix','notes');


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% TE images now...

TE_imnames = get_TE_list( image_source_folder );
mymkdir( output_mat_file_name );

fprintf('\nCopying TE files [semi-silent mode]...\n\n');

cp_cnt  = 0;
cp_step = 500;

nn = length( TE_imnames );

% make 30 images at first
%nn = 30; 

for ii = 1:nn,
 cur = TE_imnames{ii};
 
 if exist( cur,'file' ),
 % copyfile( cur,output_mat_file_name );
   load( cur );
   
   % Attention: very few of 'dummyDDB' could be empty -- avoid them!
   if isempty( dummyDDB ),continue;end
   signature_matrix = [dummyDDB; 0];
   [ju,na,ex] = fileparts( cur );
   current_savename = fullfile( output_mat_file_name,[na ex] );
   save( current_savename,'category_names','dataset_name','image_paths','signature_labels','signature_matrix','notes' );
   cp_cnt = 1 + cp_cnt;
   if cp_cnt == cp_step,
    fprintf('%04u/%u copied...\n',ii,nn);
    cp_cnt = 0;
   end % if
 end % if exist
end % ii

fprintf('.. finished ..\n');
end % eofun


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function mymkdir( name )
if exist(name,'dir'),return;end
stat = mkdir( name );
if ~stat,error(['>>mymkdir<< **cannot make folder ' name ',abort**']);end
end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function TE_imnames = get_TE_list(image_source_folder),
buffer = dir(fullfile(image_source_folder,'*.mat'));
Nbuf = length(buffer); cnt = 0; 
TE_imnames = {};
for ibuf = 1:Nbuf,
 Dname = buffer(ibuf).name; dirFlag = buffer(ibuf).isdir;
%fprintf('## %02u /%02u  %s\t%d  ',ibuf,Nbuf,Dname,dirFlag);
 if ~dirFlag, cnt = cnt + 1; TE_imnames{cnt} = fullfile(image_source_folder,Dname); end % if ~dirFlag
end % ibuf
end % eofunc



% This function goes over the image source folder 
% (which contains subfolders with original image files).
% The goal is to count those subfolders ( aka classes )
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [imgFold_Names, Ncla] = get_class_list(image_source_folder),
%%%curdir = pwd;cd(image_source_folder); unix('ls'); cd(curdir);
buffer = dir(image_source_folder);
Nbuf = length(buffer); cnt_fld = 0; imgFold_Names = {};
%%disp(Nbuf);
for ibuf = 1:Nbuf,
 Dname = buffer(ibuf).name; dirFlag = buffer(ibuf).isdir;
%fprintf('## %02u /%02u  %s\t%d  ',ibuf,Nbuf,Dname,dirFlag);
 if dirFlag & ( ~strcmp(Dname,'.') & ~strcmp(Dname,'..') ), 
  cnt_fld = cnt_fld + 1; imgFold_Names{cnt_fld} = Dname;
 %fprintf('\t\t%u  %s\n',cnt_fld,imgFold_Names{cnt_fld});
 else,  %fprintf('\n');
  continue;
 end
end
Ncla = cnt_fld;
for cla = 1:Ncla,fprintf('%s \n',imgFold_Names{cla});end
end % eofunc



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function dummyDDB2 = convert_to_one_column(dummyDDB),
[m,n] = size(dummyDDB);
if m~=3,dummyDDB2 =  dummyDDB; return; end
dummyDDB2 = []; for jj = 1:n, dummyDDB2 = [dummyDDB2 dummyDDB(:,jj)]; end % for
dummyDDB2 = dummyDDB2(:);
end % eofun


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function   siglab2 = check_labels(siglabs),
if size(siglabs,1)==3, 
 siglab2 = [siglabs(1,:) siglabs(2,:) siglabs(3,:)];
else,siglab2 = siglabs;return;end
if prod(size(siglabs)) ~= prod(size(siglab2)),error('<check_labels>: labeling is wrong, abort');end
end % eofun


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [DDB2,class_vec2,trainImgs2] = tossDDB(DDB,class_vec,trainImgs),
uc = unique(class_vec);
DDB2 = []; class_vec2 = []; trainImgs2 = []; %cell(size(trainImgs));
for ii = 1:length(uc),
 cur_list = find(class_vec==ii);
 class_vec2 = [class_vec2 class_vec(cur_list)];
 DDB2 = [DDB2 DDB(:,cur_list)];
 trainImgs2 = [trainImgs2 trainImgs(cur_list)];
end % for
end % eofun


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% read in list of classes along with number of images/tiles for each class
function [slide_names,classes] = read_info_file(inf_fname),
fid = fopen(inf_fname, 'rt'); 
tline = []; classCnt = 0;
% Reading 'reassemble_inf' file...
while feof(fid) == 0, tline = fgetl(fid); 
 if strcmp(lower(tline),'::eof'),break;end
 if strcmp(lower(tline(1:7)),':class:'),
  className = tline(8:end); classCnt = classCnt + 1; slideCnt = 0;
  classes{classCnt} = className;
 else,
  slideCnt = slideCnt + 1;
  [slide_names{classCnt}{slideCnt}, rem0] = strtok(tline);
  slide_tiles{classCnt}{slideCnt} = str2num(rem0);
 end
end % while
fclose(fid);
end % eofun

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function y = isnum(x),
y = ~isnan(str2double(x));
end % eofun


