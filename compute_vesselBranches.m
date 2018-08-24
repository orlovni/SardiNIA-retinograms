% mask_path = []; skel_path = []; VTlist = compute_vesselBranches( mask_path,skel_path );
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function VTlist = compute_vesselBranches( mask_path,skel_path ),

if isempty(mask_path),
 mask_path = fullfile('../_extracted_skeletonized','0001_test_segm_lmse.png');
end % if isempty
if isempty(skel_path),
skel_path = fullfile('../_extracted_skeletonized','0001_test_skel_lmse.png');
end % if isempty

mask = import_bw( mask_path );
skel = import_bw( skel_path );

if isempty( find(skel) ),VTlist = []; return; end


DumpYN = logical(0);
if DumpYN,dump_nm = 'dump00.txt';end % DumpYN

kernel_size = 9;
%kernel_size = 13;
%kernel_size = 17;

% maximal allowed [W]idth (eg thickness) constant:
%Wmax = 2*kernel_size-1;
Wmax = 5;
Wmax = round(.5*(kernel_size+1));
if Wmax < 7,Wmax = 7;end

% Starting point: get it from 'update_list_Ypoints'
[xstart,ystart,Yx,Yy] = update_list_Ypoints( skel );

if isempty(Yx),VTlist = []; return; end

% VesselBranchList: the vessels are segmented by the branches (from End points)
% Start the VesselBranchList:
VesselBranchList{1} = [xstart ystart];

% initiating points:
x = xstart; y = ystart;
if DumpYN,fid_dump = dump_xy(x,y,dump_nm,logical(1));end % DumpYN

keep_up = logical(1);
ifPlot  = logical(0);
orig_skel = skel; % skel is subject to degrade...
hi = [];
while keep_up,
 [x,y,keep_up,skel,hi,VesselBranchList] = crawl_vessel(skel,kernel_size,x,y,ifPlot,hi,Yx,Yy,VesselBranchList);
 if DumpYN,dump_xy(x,y,dump_nm,logical(0),fid_dump);end % DumpYN
end % while
% Remove repetitive branches if any:
VesselBranchList = cleanUp_VesselBranches( VesselBranchList );

if DumpYN,fclose(fid_dump);end % DumpYN

% recoved skel:
skel = orig_skel; clear orig_skel

% Here I run the kernel on 'unordered' mask/skeleton:
vw = run_kernel( mask,skel,kernel_size,Wmax );
VTlist = get_segmentedVThick( vw,skel,VesselBranchList );
VTlist = check_allBranchThicknesses( VTlist );
end % end function


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [x1,y1] = findAnyOpenBifurcationPoint( skel,Yx,Yy ),
%
% Given list {Yx,Yy} detects in random a point having non-empty
% neighborhood
%
x1 = []; y1 = [];

if isempty(Yx),
 error('** Problem[findAnyOpenBifurcationPoint]: empty list of Bifurc.points,abort **');
end % abort
FMIN = sqrt(2);
[y,x] = find( skel );

for ii = 1:length(Yx),
 currX = Yx(ii);
 currY = Yy(ii);
 rr = sqrt((x-currX).^2+(y-currY).^2); 
 [fmin,kmin] = min(rr);
 if fmin > FMIN,continue;end
 
 x1 = Yx(ii); y1 = Yy(ii); 
 break;
 
end % ii
end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function vessel_width = equation_constraint(t1,t2,x0,y0,xx,yy,mm,nn,Wmax,current_mask),
rr = sqrt((xx-x0).^2 + (yy-y0).^2);
ind = find( rr <= Wmax );
xx = xx(ind);
yy = yy(ind);

xx = unique(xx);
ye1 = (y0-t1*x0) + t1*xx;
ye2 = (y0-t2*x0) + t2*xx;

ind1 = find( ismember(round(ye1),1:mm) & ismember(xx,1:nn) );
ind2 = find( ismember(round(ye2),1:mm) & ismember(xx,1:nn) );

x1 = xx(ind1); y1 = ye1(ind1);
x2 = xx(ind2); y2 = ye2(ind2);

% select x1/y1: closest to the current_mask
[ym,xm] = find(current_mask);
included = [];

% pick longest of the two:
xL = x1; yL = y1;
if length(x2)>length(x1),
 xL = x2; yL = y2;
end % longest

if isempty(xL) | length(xL)<2,
 vessel_width = [];
 return;
end
for ii = 1:length(xL),
 maskdist = sqrt((xL(ii)-xm).^2 + (yL(ii)-ym).^2);
 mins(ii) = min(maskdist);
end % x1
included = find(round(mins) <= 1);
if isempty(included) | length(included)<2,
 vessel_width = [];
 return;
end
xx = xL(included);
yy = yL(included);
vessel_width = norm([xx(1);yy(1)]-[xx(end);yy(end)]);
end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function vessel_width = process_current_block( y00,x00,current_mask,current_skel,Wmax ),
vessel_width = [];

% Based on closest distance from the center of vessel to non-vessel point,
% see in foo.m

% Close (Wmax wise) to the frame center and True (mask wise):
[ny,nx] = size(current_mask);
Cx = nx/2; Cy = ny/2;

%[X,Y] = meshgrid(1:nx,1:ny); X = X - Cx; Y = Y - Cy;[F,R] = cart2pol(X,Y);[yy,xx] = find( current_skel & (R <= Wmax) );clear X Y R F;
%x0 = mean(xx); y0 = mean(yy);

[yske,xske] = find( current_skel );
R = (xske-Cx).^2 + (yske-Cy).^2; [R,ii] = sort( R );
xske = xske(ii); yske = yske(ii); x0 = xske(1); y0 = yske(1);
clear xske yske R;

if ~current_mask(y0,x0),
 [x0,y0] = get_center_fromMask( current_mask,Wmax );
end

[yk,xk] = find( ~current_mask );
% whole frame is taken by mask: width is the frame size
if isempty(yk),
 noyk = 1;
 vessel_width = size(current_mask,1);
 return;
end % isempty(yk)

rk = (xk-x0).^2 + (yk-y0).^2;
[rk,ik] = sort( rk,'ascend' );xk = xk(ik); yk = yk(ik);

% I need no more than 10 points:
n = 25;
xk = xk(1:min(length(xk),n)); yk = yk(1:min(length(xk),n)); rk = rk(1:min(length(xk),n));

% point 1: closest to the center (eg min to rk)
k0 = ik(1); xk0 = xk(1); yk0 = yk(1);

% vectors for center (c0) and the first point
c0 = [x0; y0]; v0 = [xk0; yk0];

% point 2: second closest to center && (diffX>1 OR diffY>1)
diffX = xk - xk0; diffY = yk - yk0;
cnt = 0;
for ii = 2:length(xk),
 % compute distance b/w current canditate and 'v0'
 v1 = [xk(ii);yk(ii)]; rn = norm(v1-v0);
 % compute distance b/w current canditate and 'center'
 rc = norm(v1-c0);
 if (rn>rc) & ( abs(diffX(ii))>1 | abs(diffY(ii))>1 ),
    cnt = ii;break;
 end % if
end % ii
if cnt < 2
 aaa = 1;
 vessel_width = 0; 
 return;
end
xk1 = xk(cnt); yk1 = yk(cnt);
v1 = [xk1; yk1]; rn = norm(v1-v0);

% The distance b/w the mask (white) points is lesser by 1 than the one I
% computed, the distance b/w the closest non-mask points.

vessel_width = rn - 1; 
return;


% Dead code below...
interactivePlot = logical(0);
if interactivePlot,
 figure,imagesc( current_mask ),colormap gray; set( gca,'ydir','normal' );
end % interactivePlot

[yy,xx] = find( current_skel );
x0 = mean(xx); y0 = mean(yy);

if abs( xx(1)-xx(end) ) < 1e-8,gr = 0; return;
elseif abs( yy(1)-yy(end) ) < 1e-8, gr = pi/2; return;
else,
 gr = atan( (yy(1)-yy(end)) ./ (xx(1)-xx(end)) );
end % xx = const

tg1 = tan(pi/2+atan( (yy(1)-y0) ./ (xx(1)-x0) ));
tg2 = tan(pi/2+atan( (y0-yy(end)) ./ (x0-xx(end)) ));

if interactivePlot,
 % Compute gradient ('gr'):
 if abs( xx(1)-xx(end) ) < 1e-8,gr = 0; 
 elseif abs( yy(1)-yy(end) ) < 1e-8, gr = pi/2; 
 else, gr = atan( (yy(1)-yy(end)) ./ (xx(1)-xx(end)) );
 end % Compute gr
 
 theta  = gr;
 radius = norm([x0 y0]);
 radius = 2;
 [u,v] = pol2cart( theta,radius );
 hold on,plot( [x0 x0+radius*cos(theta)],[y0 y0+radius*sin(theta)],'m' );
end % interactivePlot
clear xx yy;

[yy,xx] = find( current_mask );
[mm,nn] = size(current_mask);
vessel_width = equation_constraint(tg1,tg2,x0,y0,xx,yy,mm,nn,Wmax,current_mask);
end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [x0,y0] = get_center_fromMask( current_mask,Wmax ),
[ny,nx] = size(current_mask);
Cx = nx/2; Cy = ny/2;
[X,Y] = meshgrid(1:nx,1:ny); 
X = X - Cx; Y = Y - Cy;
[F,R] = cart2pol(X,Y);
[yy,xx] = find( current_mask & (R <= Wmax) );
x0 = mean(xx); y0 = mean(yy);
end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function VesselBranchList1 = cleanUp_VesselBranches( VesselBranchList ),
% Some of VesselBranchList come out as existing repetition of other
% branches. Here I cut those out
cnt = 0;
wholeList = 1:length( VesselBranchList );
for ii = wholeList,
 cur = VesselBranchList{ii};
 if isempty( cur ),continue;end
 
% if cur belongs to any of VesselBranchList{jj},jj ~= ii, skip ii step
 skip = logical(0);
 if size(cur,1) == 1,
  checkList = setdiff( wholeList,ii );
  x = cur(1); y = cur(2); 
  for jj = checkList,
   if ismember(x,VesselBranchList{jj}(:,1)) & ismember(y,VesselBranchList{jj}(:,2)), skip = logical(1);break;end % skip
  end % jj
 end % size(cur,1) == 1
 if skip,continue;end
 
 cnt = cnt + 1;
 VesselBranchList1{cnt} = VesselBranchList{ii};
end % ii
end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function bw = import_bw( impath ),
a  = imread( impath );
bw = logical(a);
end % eofunc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [xL,yL,xx,yy] = update_list_Ypoints( skel ),
xL = []; yL = []; xx = []; yy = [];

cn = get_all_cross_numbers( skel );

% 1( E ): ridge (E)nding point     
[y1,x1] = find( cn==1 );

% 3( Y ): bifurcation point      
[y3,x3] = find( cn==3 );

% 4( X ): crossing point         
[y4,x4] = find( cn==4 );

[yy,xx] = find( cn==3 | cn==4 );

% Priority: #1 (X), #2 (Y), #3 (E)
if ~isempty(x4),xL = x4(1);yL = y4(1); return; end
if ~isempty(x3),xL = x3(1);yL = y3(1); return; end

%if isempty(xx),error('** there''s neither crosses nor bifurcations,abort **');end
if isempty(xx),
 [yy,xx] = find( cn==1 );
 xL = xx(1);yL = yy(1);
end
return;

if ~isempty(x1),xL = x1(1);yL = y1(1); return; end

[y,x] = find( skel );
if isempty(x), return; end
end % eofunc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [limx,limy] = axLims(skel),
[y,x] = find( skel );
limx = [min(x) max(x)];
limy = [min(y) max(y)];
end % eofunc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [x,y,keep_up,skel,hi1,VesselBranchList] = crawl_vessel(skel0,kernel_size,x0,y0,ifPlot,hi0,Yx,Yy,VesselBranchList0)
%
% Dynamically updates list of coordinates (x,y) of the vessel ridge. 
% Once the list is empty, the ridge run stops.
%
keep_up = logical(1);
skel    = skel0;
hi1 = hi0;
VesselBranchList = VesselBranchList0; 
clear VesselBranchList0;

[ma,na] = size( skel );

if ifPlot & isempty(hi1) & ~isempty(x0),
 figure; hi0 = imagesc(skel);colormap gray;set( gca,'ydir','normal' );
 hold on; plot(Yx,Yy,'om','markersize',16,'linewidth',2);
 %[limx,limy] = axLims(orig_skel); xlim(limx); ylim(limy);
end

[x,y,skel,contX,hi1,VesselBranchList] = ifCloseToBorder( kernel_size,x0,y0,skel,ifPlot,hi1,Yx,Yy,VesselBranchList );
if contX,return;end
  
 % if normal point: delete it, find the next closest to it
skel(y,x) = logical(0);
if ifPlot,
  [x,y,hi1,VesselBranchList] = closing_dance(x,y,skel,ifPlot,hi1,Yx,Yy,VesselBranchList);
else,
  [x,y,foo,VesselBranchList] = closing_dance(x,y,skel,ifPlot,[],Yx,Yy,VesselBranchList);
end % ifPlot
 
 % update the stop sign
[yy,xx] = find(skel);nx = length(xx);
if isempty(x) | isempty(x0) | nx < 2,
  keep_up = logical(0);
  if ifPlot,close(gcf);end
end % update the stop sign
end % eofunc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [x1,y1,skel2,contX,varargout] = ifCloseToBorder( kernel_size,x,y,skel,ifPlot,varargin ),
skel2 = skel; x1 = x; y1 = y;
[ma,na] = size( skel ); 
varargout{1} = [];
varargout{2} = [];
contX = logical(0); % if true, then 'continue' on X (or Y) is execulted

Yx = varargin{2};
Yy = varargin{3};
VesselBranchList = varargin{4};

if ifPlot,
 hi = varargin{1}; hi1 = hi; varargout{1} = hi1;
else, hi = [];
end % ifPlot

if x-kernel_size < 2 | x+kernel_size > na-1,
  contX = logical(1);
  skel2(y,x) = logical(0);
  [x1,y1,hi1,VesselBranchList] = closing_dance(x,y,skel,ifPlot,hi,Yx,Yy,VesselBranchList);
  if ~isempty(x1),VesselBranchList{end} = [VesselBranchList{end}; x1 y1];end
  if ifPlot,
   varargout{1} = hi1;
  end % ifPlot
end % if x-kernel_size

if y-kernel_size < 2 | y+kernel_size > ma-1,
  contX = logical(1);
  skel2(y,x) = logical(0);
  [x1,y1,hi1,VesselBranchList] = closing_dance(x,y,skel,ifPlot,hi,Yx,Yy,VesselBranchList);
  if ~isempty(x1),VesselBranchList{end} = [VesselBranchList{end}; x1 y1];end
  if ifPlot,
   varargout{1} = hi1;
  end % ifPlot
end % if x-kernel_size
varargout{2} = VesselBranchList;
end % eofunc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function VTlist = get_segmentedVThick( vw,skel,VBlist ),
%
% Computes width of all vessels by running the kernel of 'kernel_size' width
%
[ma,na] = size(skel);
% detect skeleton whites:
[ys,xs] = find( skel );
xys = [xs(:) ys(:)];

% Given xy:
xyg = vw(:,1:2);
VTlist = [];
Nbranches = length(VBlist);

for k = 1:Nbranches,
    
 nk = length(VBlist{k});

% current branch:
 xy = VBlist{k}(:,1:2);
 [foo,jj] = ismember(xy,xyg,'rows'); %aa = xyg(:); foo2 = [aa(jj(:,1)) aa(jj(:,2))]; clear foo aa
 foo2 = xyg(jj,:);
 criterion = sum(sum(foo2-xy))==0; clear foo2 xy

% if criterion, then xy is the branch:
 if criterion,
  VTlist{k} = vw(jj,:);
 else,
  VTlist{k} = ones(1,3);
 end

end % k, 1:Nbranches
end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function vw = run_kernel( mask,skel,kernel_size,Wmax ),
%
% Computes width of all vessels by running the kernel of 'kernel_size' width
%
[ma,na] = size(mask);
% detect skeleton whites:
[iarray,jarray] = find( skel );
n = length(iarray);
vw = [];
for ii = 1:n,
 y = iarray(ii);
 x = jarray(ii);
 if x-kernel_size < 2 | x+kernel_size > na-1,vw = [vw; x y 0];continue;end
 if y-kernel_size < 2 | y+kernel_size > ma-1,vw = [vw; x y 0];continue;end
 
 % process current block
 xblock = x-kernel_size:x+kernel_size;
 yblock = y-kernel_size:y+kernel_size;
 vessel_width = process_current_block( y,x,mask(yblock,xblock),skel(yblock,xblock),Wmax );
 if isempty(vessel_width),vw = [vw; x y 0];continue;end
 vw = [vw; x y vessel_width];
 
end % ii

end % eofunc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function VT = check_allBranchThicknesses( VTlist ),
VT = VTlist;
n = length(VTlist);
for ii = 1:n,
 t = VTlist{ii}(:,3);
 k = find(t==0);
 if ~isempty(k),
  try
   t = interpolate_zero_values( t );
  catch,
  %whos t
   t(t==0) = 1;
  end % try
 %VT{ii} = [VTlist{ii}(:,1:2) t];
  VT{ii}(:,3) = t;
 end % if ~empty
end % ii
end % eofunc


% interpolate zero values with the closest neighbor
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function t2 = interpolate_zero_values( t )
n = length(t);
x = (1:n)';
x0 = find(t==0);
y0 = zeros(size(x0));
x1 = find(t>0);
y1 = t(x1);

if length(x1)==1 & length(x0)>1,
 t2 = mean(t).*ones(size(t)); 
 return;
end % if
y0 = interp1(x1,y1,x0,'spline');
y0(y0<0) = 1;
t2 = zeros(size(t)); t2(x0) = y0; t2(x1) = y1;
end % eofunc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [x1,y1,hi1,VesselBranchList] = closing_dance(x,y,skel,ifPlot,varargin),
%
% Restarts the ridge run.
% Ends/concludes ongoing ridge run and then restarts from (x1,y1).
% The point (x1,y1) is found with 'findAnyOpenBifurcationPoint'.
%
FMIN = sqrt(2);
hi = [];
%varargout{1} = [];
hi1 = [];
Yx = varargin{2}; 
Yy = varargin{3}; 
VesselBranchList = varargin{4};

if ifPlot, 
 hi = varargin{1}; 
%orig_skel = varargin{4};
end % ifPlot

% Find the next point (x1,y1) as the one closest to (x,y)
[cy,cx] = find( skel ); 

% if x1==x & y1==y this is wrong: exclude this point from the lists cx,cy
if ~isempty(cx),
 rr = sqrt((x-cx).^2+(y-cy).^2); [fmin,kmin] = min(rr);
 if fmin==0,allIND = 1:length(rr);
  kleans = setdiff( allIND,kmin );
  cx = cx(kleans); cy = cy(kleans);
 end % fmin==0
end % ~isempty(cx)

if ~isempty(cx), 
    
rr = sqrt((x-cx).^2+(y-cy).^2); 
[fmin,kmin] = min(rr);
x1 = cx(kmin); y1 = cy(kmin);

if fmin <= FMIN,
 if ifPlot,
  delete(hi);hi1 = imagesc(skel);colormap gray;set( gca,'ydir','normal' );
  hold on; plot(Yx,Yy,'om','markersize',16,'linewidth',2);
  %[limx,limy] = axLims(orig_skel); xlim(limx); ylim(limy);
  drawnow;
  %varargout{1} = hi1;
 end % ifPlot
 VesselBranchList{end} = [VesselBranchList{end}; x1 y1];
 return;
end % fmin <= FMIN

else, 
 %x1 = []; y1 = []; if ifPlot,varargout{1} = [];end
 %return;
end % ~isempty(cx)

% If fmin > 1 (eg if distance is larger than normally) then the next start
% point is not the "closest to the current x,y" (that is FAR) but instead
% any Y point. If there's no Y points, start from an end point
% Thus, I am looking for Y/end points here.
% However, the concept seems ill-conceived: update would erase the
% bifurcation (Y/end) points as well.
%%[x1,y1] = update_list_Ypoints( skel );
%
% A better concept is to pick Any of bifurcation (Y/end) points and check
% whether its vicinity contains remaining ridge points to run. If it does,
% run from there. Otherwise, try the next bifurcation point; etc.
x1 = Yx; y1 = Yy;
[x1,y1] = findAnyOpenBifurcationPoint( skel,Yx,Yy );

% Update the VesselBranchList:
VesselBranchList{end+1} = [x1 y1];
%varargout{2} = VesselBranchList;

if ifPlot,
%while ~strcmpi( inputdlg('Next Y point? (y/n)','Continue?',1,{'y'}),'y' ),continue; end % while 
 delete(hi);hi1 = imagesc(skel);colormap gray;set( gca,'ydir','normal' );
 hold on; plot(Yx,Yy,'om','markersize',16,'linewidth',2); drawnow;
%varargout{1} = hi1;
end % ifPlot
end % eofunc


%
% 'get_all_cross_numbers.m' is another core part of
% 'foo_crossing_number.m' to be used from within other projects.
%
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
function yn = isodd( x ),r = rem( x,2 );
if r > 0,yn = logical(1);
elseif r==0,yn = logical(0);
end
end % eofunc
function yn = iseven( x ),yn = logical(1);
if isodd( x ),yn = logical(0);end
end % eofunc

