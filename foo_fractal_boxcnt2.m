%
% foo_fractal_boxcnt2( '../__computechains/__DifferentImages/01_manual1-Bess.tif' );
% foo_fractal_boxcnt2( '../__computechains/__DifferentImages/01_manual1.gif' );
%
% foo_fractal_boxcnt2( '../__computechains/__swSkel/lgFFT_segm512.tif' );
% foo_fractal_boxcnt2( '../__computechains/__swMask/lgFFT_segm512.tif' );
% foo_fractal_boxcnt2( '../__computechains/__ManualSkeleton/lgFFT_orig.tif' );
% foo_fractal_boxcnt2( '../__computechains/__ManualMask/lgFFT_orig.tif' );
% foo_fractal_boxcnt2( '../__computechains/__OrigImage/lgFFT_orig.tif' );
%
%
function foo_fractal_boxcnt2( imname ),

if exist( imname,'file' ),a = imread( imname );end
%disp(size(a));
if length(size(a)) > 2,a = a(:,:,1);end

% assuming image is uint8:
a = im2double(a);
a = ones(size(a))-a; % flip foreground/background

% Gray to BW
th = 0.5;
if ~isempty( find(a>0 & a<1) ),
 a(a>th) = 1; a(a<=th) = 0;
end % ~isempty

%figure,imagesc(a),colormap gray;
d = hausDim( a );
fprintf('%s: %.6f\n',imname,d);
end % eofunc



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function y = hausDim( I )
%function [ D ] = hausDim( I,imname )
% HAUSDIM Returns the Haussdorf fractal dimension of an object represented by
% a binary image.
%    Returns the Haussdorf fractal dimension D of an object represented by the
%    binary image I. Nonzero pixels belong to an object and 0 pixels 
%    constitute the background.
%
%    Algorithm
%    ---------
%    1 - Pad the image with background pixels so that its dimensions are a 
%        power of 2.
%    2 - Set the box size 'e' to the size of the image.
%    3 - Compute N(e), which corresponds to the number of boxes of size 'e' 
%        which contains at least one object pixel.
%    4 - If e > 1 then e = e / 2 and repeat step 3.
%    5 - Compute the points log(N(e)) x log(1/e) and use the least squares 
%        method to fit a line to the points.
%    6 - The returned Haussdorf fractal dimension D is the slope of the line.
%
%    Author
%    ------
%    Alceu Ferraz Costa 
%    email: alceufc [at] icmc [dot] usp [dot] br
%

    % Pad the image with background pixels so that its dimensions are a power of 2.
    maxDim = max(size(I));
    newDimSize = 2^ceil(log2(maxDim));
    rowPad = newDimSize - size(I, 1);
    colPad = newDimSize - size(I, 2);
    
%   I = padarray(I, [rowPad, colPad], 'post');
%   I = array_padd(I, [rowPad, colPad],0, 'post');
    if 1==2
    I = [zeros(size(I,1),colPad) I zeros(size(I,1),colPad)];
    I = [zeros(rowPad,size(I,2)); I; zeros(rowPad,size(I,2))];
    end % 1==2

    rPadUp = floor(rowPad/2); rPadDo = rPadUp; 
    if sum([rPadUp rPadDo size(I, 1)]) < newDimSize,rPadDo = rPadDo+1;end
    
    cPadLe = floor(colPad/2); cPadRi = cPadLe; 
    if sum([cPadLe cPadRi size(I, 2)]) < newDimSize,cPadRi = cPadRi+1;end
    
    I = [zeros(size(I,1),cPadLe) I zeros(size(I,1),cPadRi)];
    I = [zeros(rPadUp,size(I,2)); I; zeros(rPadDo,size(I,2))];
   

    boxCounts = zeros(1, ceil(log2(maxDim)));
    resolutions = zeros(1, ceil(log2(maxDim)));
    
    boxSize = size(I, 1);
    boxesPerDim = 1;
    idx = 0;
    while boxSize >= 1
        boxCount = 0;
        
        for boxRow = 1:boxesPerDim
            for boxCol = 1:boxesPerDim
                minRow = (boxRow - 1) * boxSize + 1;
                maxRow = boxRow * boxSize;
                minCol = (boxCol - 1) * boxSize + 1;
                maxCol = boxCol * boxSize;
                
                objFound = false;
                for row = minRow:maxRow
                    for col = minCol:maxCol
                        if I(row, col)
                            boxCount = boxCount + 1;
                            objFound = true; % Break from nested loop.
                        end;
                        
                        if objFound
                            break; % Break from nested loop.
                        end;
                    end;
                    
                    if objFound
                        break; % Break from nested loop.
                    end;
                end;
            end;
        end;
        
        idx = idx + 1;
        boxCounts(idx) = boxCount;
        resolutions(idx) = 1 / boxSize;
        
        boxesPerDim = boxesPerDim * 2;
        boxSize = boxSize / 2;
    end;
    
    %D = polyfit(log(resolutions), log(boxCounts), 1); D = D(1);
    
    %slopes = piecewise_polyfit( log(resolutions),log(boxCounts),imname  );
    bn = log2(boxCounts); bn = boxCounts./sum(boxCounts);
    
    % Vector of Housdorff dimensions by scales (N.O. Mar 2016):
    x  = resolutions./resolutions(1);
    x( x==1 ) = 1 + eps;
    y  = log(boxCounts)./log(x);
   %figure,hs=loglog(x(2:end),y(2:end),'r-');grid on;set(hs,'linewidth',1.5);xlim([1 2e3]);xlabel('x');ylabel('log(N)/log(x)');
    y = y(2:end);
return;    
%foo = exp(slopes.^2); foo = foo - min(foo);
fprintf('%s\n',num2str(bn,' %.5e'));
fprintf('%s\n',num2str(y ,' %.5e'));
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function slopes = piecewise_polyfit( x,y,varargin ),

if length(varargin) > 0,
 foo = varargin{1};
 [jnk,ttle,ex] = fileparts(foo);
 ifprint = logical(1);
else,ttle = []; ifprint = logical(0);
end % if title exists

n = length(x)-1; 
%n = 6;
slopes = zeros(1,n);
for ii = 1:n,
 x0 = x(ii); x1 = x(ii+1); dx = x1-x0;
 y0 = y(ii); y1 = y(ii+1); dy = y1-y0;
 if abs(dx) <= eps, slopes(ii) = sign(dy)*sign(dx);
 elseif abs(dx) > eps, slopes(ii) = dy / dx;
 end % if
end % ii
if ~ifprint,return;end

%disp(mean(slopes(1:n)));
fprintf('%s\n',num2str(mean(slopes(1:n)),' %.4f'));
%fprintf('%s\n',num2str(slopes,' %.3f'));
%fprintf('%s\n',num2str(log2(slopes),' %.3f'));

%foo = (10+slopes).^2; foo = foo - min(foo);
foo = exp(slopes.^2); foo = foo - min(foo);
fprintf('%s\n',num2str(foo,' %.2f'));
return;

figure,plot(x,y,'-ob');hold on;axis([min(x) max(x) min(y) max(y)]);
plot([x(1) x(end)],[y(1) y(end)],'-r');title(ttle,'Interpreter','none');
end % eofunc
