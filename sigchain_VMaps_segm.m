%
% imSG = []; imSk = []; [signature_vector,sigSizes] = sigchain_VMaps_segm( imSG,imSk );
%
% Chain for segmented vessel maps includes features:
% Fractal dimension, Thickness histogram
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Nikita Orlov, norlov@nih.gov
% Aug 2016 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
function [signature_vector,sigSizes,isaborted] = sigchain_VMaps_segm( mask_path,skel_path,masterCaliber ),

if isempty( mask_path ),
 mask_path = fullfile('../_extracted_skeletonized','0001_test_segm_lmse.png');
end % if isempty
if isempty( skel_path ),
 skel_path = fullfile('../_extracted_skeletonized','0001_test_skel_lmse.png');
end % if isempty

isaborted = logical(0);

% generate list of vessicles containing [x(:) y(:) thickness(:)]...
XYT_by_branches = compute_vesselBranches( mask_path,skel_path );

% no ridge points/crosses/bifurcations:
if isempty(XYT_by_branches),
 isaborted = logical(1);
 signature_vector = []; 
 sigSizes = [];
 record_mean_caliber( mask_path,[],isaborted );
 return;
end

%save wrk_branches XYT_by_branches

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% recast cell array 'XYT_by_branches' into a struc array 'branch()'
%load wrk_branches;

% evaluate thickness values for all branches in given map
allThick = evaluate_thickness( XYT_by_branches );


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Nota Benne:
% 'uberMean' is an important measure of overall boderline between thinner
% and thicker vessels. It is based on a representative number of different
% vessel maps.
%uberMean = mean( allThick );
uberMean = masterCaliber;

% If the mean caliber does not exist: save it and break out w/o computing
% sig vector. If the mask mean caliber is detected: it must be read earlier
if isempty( uberMean ),
 record_mean_caliber( mask_path,mean( allThick ),isaborted );
 signature_vector = [];
 sigSizes = [];
 return;
end % isempty

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
branch = initiate_branch( XYT_by_branches,uberMean );

% compute whole-map Tortuosity (15-vec)
tortuosity_WholeMap = get_tortuosity( branch );
signature_vector = tortuosity_WholeMap;
sigSizes = length( tortuosity_WholeMap );
L15 = 15;

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% compute whole-map Fractal Dimensions, 8-vec
[FD_wolemap,xWhMap,yWhMap] = get_FractalDim( branch,logical(1) );
signature_vector = [signature_vector FD_wolemap];
sigSizes = [sigSizes length(FD_wolemap)];
L8 = 8;

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% compute whole-map Junction analysis, 3-vec
JA_wholemap = get_JunctionAnal( xWhMap,yWhMap,uberMean );
signature_vector = [signature_vector JA_wholemap];
sigSizes = [sigSizes length(JA_wholemap)];
L3 = 3;

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Segregate by thickness. Each Thickness Category is a separate channel
brCategory = THick_stratify( branch );

nCh = length( brCategory );
for ch = 1:nCh,
 currCategory = brCategory{ch};
 
 % Compute tortuosity (15-vec) for given thickness category
 if isempty( currCategory ),
  tortuosity{ch} = zeros(1,L15);
 else,
  tortuosity{ch} = get_tortuosity( currCategory );
 end % isempty(currCategory)
%sigSizes = [sigSizes length(tortuosity{ch})];
 sigSizes = [sigSizes L15];
 signature_vector = [signature_vector tortuosity{ch}];
 
 % Compute branch-average Fractal Dimensions (8-vec) for given category
 if isempty( currCategory ),
  FD_brAvg{ch} = zeros(1,L8);
 else,
  FD_brAvg{ch} = get_FractalDim( currCategory,logical(0) );
 end % isempty(currCategory)
%sigSizes = [sigSizes length(FD_brAvg{ch})];
 sigSizes = [sigSizes L8];
 signature_vector = [signature_vector FD_brAvg{ch}];
 
 % Compute Junction analysis (3-vec) for given thickness category
 if isempty( currCategory ),
  JA{ch} = zeros(1,L3);
 else,
  [xx,yy] = get_xy_fromBranch( currCategory );
  JA{ch} = get_JunctionAnal( xx,yy,uberMean );
 end % isempty(currCategory)
%sigSizes = [sigSizes length(JA{ch})];
 sigSizes = [sigSizes L3];
 signature_vector = [signature_vector JA{ch}];
end % ch

end % end function


function [xx,yy] = get_xy_fromBranch( branch ),
xx = []; yy = [];
n = length( branch );
for ii = 1:n,
 xc = branch(ii).x; yc = branch(ii).y;
 xx = [xx xc']; 
 yy = [yy yc'];
end % ii
end % end function


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function record_mean_caliber( mask_path,avgThick,isaborted ),
meanCaliberFname = gen_AVGthick_fname( mask_path,isaborted );
save( meanCaliberFname,'avgThick' );
return;

SG = '_segm_';
[fileDir,fileNm,fileExt] = fileparts( mask_path );
k0 = strfind( fileNm,SG ) - 1;
k1 = k0 + length(SG);
part1 = fileNm(1:k0);
part2 = fileNm(k1:end);
meanCaliberFname = [part1 '_AVGthick_' part2];
meanCaliberFname = fullfile( fileDir,meanCaliberFname );
save( meanCaliberFname,'avgThick' );
end % end function


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


% compute whole-map Junction analysis, 3-vec
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function JA_wholemap = get_JunctionAnal( xWhMap,yWhMap,uberMean ),

% form 'skel':
skel = zeros((max(yWhMap)-min(yWhMap)),(max(xWhMap)-min(xWhMap)));
nx = length(xWhMap);
for ii = 1:nx, skel(yWhMap(ii),xWhMap(ii)) = 1; end % ii

skel = logical(skel);
cn = get_all_cross_numbers( skel );

% 1( E ): ridge (E)nding point     
[y1,x1] = find( cn==1 );

% 3( Y ): bifurcation point      
[y3,x3] = find( cn==3 );

% 4( X ): crossing point         
[y4,x4] = find( cn==4 );


% 'too close': 
% Number of crossing points is 'unstable', as it relies on skeletonization
% of the segmentation map. Skeletonization algoritm is known for its
% artificial 'bridges', when instead of a single (and true) crossing point
% it produces instead two (but close) bifurcation (Y) points. There is no
% exact algorithm to treat this defect. Given the max thickness of given
% vessel segmented map, one can say that two Y points closer than the max
% thickness should be combined together into single X point.

too_close = uberMean;
if length(x3) > 2,
 [x3,y3,x4,y4] = inspect_Ypoints( x3,y3,x4,y4,too_close );
end

JA_wholemap = [length(y1) length(y3) length(y4)];
end % end function

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [X3,Y3,X4,Y4] = inspect_Ypoints( x3,y3,x4,y4,too_close )
% x3,y3: bifurcation points ( Y )
% x4,y4: crossing points ( X )
% Compute pairwise distances b/w all bifurcation points
n = length(x3); 
X3 = []; Y3 = []; 
X4 = []; Y4 = []; 
excluded = [];
for ii = 1:n,
%   fprintf('>inspect_Ypoints< %u/%u\n',ii,n);
 curr_x = x3(ii);
 curr_y = y3(ii);
 curr_v = [curr_x; curr_y];
 restind = setdiff( 1:n,ii );
 % rest of matrix: all bifurcation points but the current. Columns
 % correspond to different bifurcation points
 restmat = [x3(restind)'; y3(restind)'];
 
 % get [minimal of course] distance from current vector to restmat
 [md,x0,y0] = get_dist2restmat( curr_v,restmat );
 
 % criterion to merge the two bifurcation points
 if md > 0 & md < too_close,
  % point 1: ii -- {x3(ii),y3(ii)}, and point 2: {x0,y0}
  i0 = find( x0==x3 & y0==y3 );
  
  % Exclude entries i0,ii from x3,y3 (if they dont exist there already)
  if isempty(intersect( excluded,[i0 ii] )),
   excluded = [excluded i0 ii];
   X4 = [X4 mean([x0 curr_x])];
   Y4 = [Y4 mean([y0 curr_y])];
  end % ~ismember
 end % condition is met
end % ii
excluded = unique( excluded );

% make exclusions
indices_remaining = setdiff( 1:n,excluded );
X3 = x3( indices_remaining );
Y3 = y3( indices_remaining );
end % eofunc


% get [minimal of course] distance from current vector to restmat
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [md,x0,y0] = get_dist2restmat( curr_v,restmat )
[m,n] = size( restmat );
vd = zeros(1,n);
for jj = 1:n, vd(jj) = norm( curr_v - restmat(:,jj) );end % jj
 [md,mj] = min( vd );
 x0 = restmat(1,mj); y0 = restmat(2,mj);
end % eofunc


% a core part of foo_crossing_number.m used from within other projects
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function cn = get_all_cross_numbers( mask ),
[im,jm] = find( mask ); 
n = length(im); ind = 1:n; cn = zeros( size(mask) );
for t = ind, ii = im(t); jj = jm(t);
 if min(ii,jj) < 2 | ii == size(mask,1) | jj == size(mask,2),continue;end
 an = get3x3neibors( double(mask),ii,jj,3 );
% get bifurcation class
 cn(ii,jj) = cross_number( mask,ii,jj ); 
end % t
end % eofunc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function an = get3x3neibors( a,ii,jj,n )
if isodd(n),n2 = (n-1)/2;else,n2=n/2;end
i3 = ii-n2:ii+n2;
j3 = jj-n2:jj+n2;
indi = 1:n; indj = 1:n;
for iy = indi,
 ix = indj;
 an(iy,ix) = a(i3(iy),j3);
end
end % eofunc
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function fs = cross_number( mask,ii,jj ),
an = get3x3neibors( double(mask),ii,jj,3 );
an = an(:)';
target_order = [8 7 4 1 2 3 6 9 8];
foo = an( target_order ); 
f1 = foo(1:end-1);
f2 = foo(2:end);
fa = abs( f1 - f2 );
fs = sum( fa )*0.5;
end % eofunc
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function yn = isodd( x ),r = rem( x,2 );
if r > 0,yn = logical(1);
elseif r==0,yn = logical(0);
end
end % eofunc
function yn = iseven( x ),yn = logical(1);
if isodd( x ),yn = logical(0);end
end % eofunc


% compute Fractal Dimension as a mean of all given branches
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [FD,x,y] = get_FractalDim( branch,isWholeMap ),
% FractalDim requires a pixel image (with all the branches given). So, I
% have to generate a tempo image composed of the branches of given
% thickness.
n = length( branch );
% maximal allowed scale is 8:
maxScale = 8;
FD_matrix = zeros(n,maxScale);
x = []; y = [];
for ii = 1:n,
%fprintf( '#%03u/%u\n',ii,n );
xx = []; yy = [];
xx = branch(ii).x'; yy = branch(ii).y';

% FD for givan branch:
if ~isWholeMap,
 bar = FractalDim_8vec( xx,yy );
 if length(bar) ~= maxScale, bar = bar(1:maxScale);end
 FD_matrix(ii,:) = bar;
end % if ~isWholeMap

x = [x xx]; y = [y yy];
end % ii

% Fractal Dimension, branch-average (8-vec):
if ~isWholeMap,
 FDmean  = mean( FD_matrix,1 ); FD = FDmean; 
 return;
 
% Fractal Dimension, whole skeleton (8-vec):
else,
 FDwhole = FractalDim_8vec( x,y );
 if length(FDwhole) ~= maxScale,FDwhole = FDwhole(1:maxScale);end
 FD = FDwhole;
 return;
end % if ~isWholeMap

end % end function


% compute Tortuosity as a mean of all given branches
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function tortuosity = get_tortuosity( branch ),
%branch = struct('x','y','thick','TortuosityLib_15vec','FractalDim_10vec');
n = length( branch );
Tortuosity_matrix = zeros(n,15);
for ii = 1:n,
%fprintf( '#%03u/%u\n',ii,n );
 foo = tortuosity_15vec( branch(ii).x,branch(ii).y );
%branch(ii).TortuosityLib_15vec = foo;
 Tortuosity_matrix(ii,:) = foo;
end % ii

% whole-map tortuosity:
tortuosity = mean( Tortuosity_matrix,1 );
end % end function



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function brCategory = THick_stratify( branch ),
THcategories = [branch.thickcategory];
uth = unique( THcategories );
for ii = unique( THcategories ),
 current = find( ii == THcategories );
 brCategory{ii} = branch(current);
end % ii
end % end func


% evaluate thickness values for all branches
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function allThick = evaluate_thickness( XYT_by_branches ),
n = length( XYT_by_branches );
for ii = 1:n,
 % whole-branch average thickness:
 allThick(ii) = mean( XYT_by_branches{ii}(:,3) );
end % ii
if mean(allThick) < 1
    aaa = 1;
end % if
end % end func

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function branch = initiate_branch( XYT_by_branches,uberMean ),
%branch = struct('x','y','thick','TortuosityLib_15vec','FractalDim_10vec');
n = length( XYT_by_branches );
for ii = 1:n,
 branch(ii).x = XYT_by_branches{ii}(:,1);
 branch(ii).y = XYT_by_branches{ii}(:,2);
 branch(ii).thick = XYT_by_branches{ii}(:,3);
 isthick = mean( XYT_by_branches{ii}(:,3) ) > uberMean/2;
 % Thickness category
 if mean( XYT_by_branches{ii}(:,3) ) <= uberMean/2, thickcategory = 1;
 else, thickcategory = 2;
 end % if
 branch(ii).thickcategory = thickcategory;
end % ii
end % end func





%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function y = norm0max( x ), x0 = min(x(:)); x1 = max(x(:));
y = x1.* (x-x0)./(x1-x0+eps);
end % eofunc
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function y = norm01( x ), x0 = min(x(:)); x1 = max(x(:));
y = (x-x0)./(x1-x0+eps);
end % eofunc



%
% helper function combines multiple function outputs into a single output vector
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function concat = concat_outputs (varargin)
N = nargin;
concat = [];
for i=1:N
	vec = double(varargin{i}); % promote to doubles
	% some vectors are row-oriented others are column-oriented.
	% This fixes them all to be row oriented
	vec = vec(:)';
	concat = [concat vec];
end
end % end function
