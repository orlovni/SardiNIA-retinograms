function [DDB,siglabs,eig_vecs,mean_img] = getSigs(imgNames,dataNameStamp)

%function [DDB,siglabs,eig_vecs,mean_img] = getSigs(imgNames,varargin)
%   INPUT:          imgNames    - filenames of input images.  i.e. 'image1.tif'
%   OUTPUT:         DDB          - matrix of signatures.  columns are images, rows are signatures.
%                                   so, DDB(1,1) is the first image's first signature and DDB(5,1) is
%                                   the first image's fifth signature.  DDB(1,5) is the fifth image's
%                                   first signature
%                   siglabs     - names of all of the signatures outputed in DDB
%                   eig_vecs    - eigen_vectors used
%                   mean_img    - mean image used
%   NOTES:
%       take in some image names, return images' signatures and then some.  eigen vectors are
%       returned because during testing, you cannot use the test images to create the eigenfaces.
%       that's also why there is an option to input the eigenfaces and the mean face.  
%   IF HACKING: 
%       please remember that 'sig_number' has been hardcoded and not dynamically assigned (running
%       out of time here...)  So, if you try and change which signatures are used, you're gonna have
%       to change that number yourself!
%   Written by Lawrence David - 2003.  lad2002@columbia.edu
%
%   Revisions by Nikita Orlov - 2003-2005.  

[imgSize imgSize] = size(imread(imgNames{1}));
eig_vecs = [];
mean_img = [];
img_quantity      = length(imgNames);         % how many images

if ~exist(dataNameStamp),
 stmkdir = mkdir(dataNameStamp);
if ~stmkdir,
 if7 = exist(fullfile(pwd,dataNameStamp));
 if if7 ~= 7, error('!! can''t create folder; folder does''t exist !!'); end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get signatures                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cnt = 0;
for i = 1:img_quantity, cnt = cnt + 1;     
    tic
    fprintf('image %03u / %03u',i,img_quantity);
    %% no lofty params:
    
    [DDBbufer,sigSizes] = typicv20(imgNames{i});
    if i == 1,
     [mddb,nddb] = size(DDBbufer); % memorize the right size of signature vec.
    else,
     if size(DDBbufer,1) ~= mddb, % if problem with sig vector...
      fprintf('  ( skipped )\n');
      cnt = cnt - 1;
      continue;
     end
    end
    DDB(:,cnt) = DDBbufer;
    
    
    timeTaken = toc;
    fprintf('  (%6.2f seconds )\n',timeTaken);
    
% dump computed file as *.mat
    if i==1, siglabs = get_siglabs(sigSizes); end
    dummyNm = imgNames{i};
    [bufDir,bufFname,bufExt] = fileparts(dummyNm); 
    dummyDDB = DDB(:,cnt);
    
    MATname = fullfile(pwd,dataNameStamp,[bufFname '.mat']);
    
    fprintf('dummyNm\t\t%s\n',dummyNm);
    fprintf('bufDir\t\t%s\n',bufDir);
    fprintf('bufFname\t\t%s\n',bufFname);
    fprintf('pwd\t\t%s\n',pwd);
    fprintf('dataNameStamp\t\t%s\n',dataNameStamp);
    fprintf('MATname\t\t%s\n',MATname);
    
    error('stop');
    
    save(MATname,'dummyDDB','siglabs','dummyNm','-V6');
end
return;


% Creating all of the labels for all of the signatures.  The idea being using "place" and "count" is
% that you don't have to modify the 'siglabs' data structure whenever you want to introduce a new
% set of signatures and more especially, when you want to remove the old one.  This way, it scales
% on its own

function siglabs = get_siglabs(sigSizes),

% revised signatures...
uniqueLabNames = {'blob size ','region sigs ','edge sigs ','hull sigs ', ...
        'Haralick end column ','Haralick mean ','Zernike ','Zernike(FFT) ', ...
        'Haralick(FFT) end column ','Haralick(FFT) mean ',... 
        'Haralick(WL) eCol ', 'Haralick(WL) mean ', 'Haralick(WLfft) eCol ', 'Haralick(WLfft) mean ',...
        'Haralick(Cheb) eCol ', 'Haralick(Cheb) mean ', 'Haralick(Cheb(FFT)) eCol ', 'Haralick(Cheb(FFT)) mean ',...
        'Cheb short ', 'Cheb(FFT) short ', 'Cheb-Forier ', 'Cheb-Forier(FFT) ', ...
        'Gabor ratio ', 'Gabor hist ', ...
        'multi-scale-hist ', 'm-sc-hist (WL) ', 'm-sc-hist (WL(FFT)) ', 'm-sc-hist (Cheb) ', 'm-sc-hist (Cheb(FFT)) ', ...
        'Radon transform ', 'Radon(FFT) ', 'Radon(WL) ', 'Radon(Cheb) ', ...
        '4-comb moments ', '4-comb moments(FFT) ', '4-comb moments(WL) ', '4-comb moments(Cheb) ', ...
        'Tamura set ', 'Tamura(FFT) ', 'Tamura(WL) ', 'Tamura(Cheb) '};

siglabs = [];
for cnt = 1:length(sigSizes),
 siglabs = addSigLabel(siglabs,uniqueLabNames{cnt},sigSizes(cnt));
end
siglabs = siglabs';
return;

function siglabs = addSigLabel(siglabs0,pattern,sigSize), siglabs = siglabs0;
buff = []; for ii = 1:sigSize, buff = [buff; [pattern sprintf('%2i',ii)]]; end
if ~isempty(buff), buff = cellstr(buff);
 else, buff = {pattern};
end
if isempty(siglabs), siglabs = buff;
else, siglabs(size(siglabs,1)+1:size(siglabs,1)+size(buff,1)) = buff;
end
return;

function print_SigLab(n1,n2,sigLab),
for ii = n1:n2, fprintf('%s; ',sigLab{ii}); end
fprintf('\n\n');

