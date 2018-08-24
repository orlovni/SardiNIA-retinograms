% Chain computed on following transforms:
% [FFT, Ch, WL, Gr, Lap, GrGr] im
% [FFT, Ch, WL, Gr, Lap] Gr
% [FFT, Ch, WL, Gr, Lap] Lap
%
% Nikita Orlov
% 01-18-2008 
%
function [signature_vector,sigSizes] = sig_chain_12transforms(im),

addpath(genpath('/iicbu/d1iicbu/orlovni/myCVScheckout/OME/Transforms'));
addpath(genpath('/iicbu/d1iicbu/orlovni/myCVScheckout/OME/Statistics'));
addpath(genpath('/iicbu/d1iicbu/orlovni/myCVScheckout/OME/Segmentation'));
addpath(genpath('/iicbu/d1iicbu/orlovni/myCVScheckout/OME/Maths'));


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Compute image transforms of the image
[FFT_imTSF,Ch_imTSF,Bess_imTSF,WL_imTSF,Grad_imTSF,Lap_imTSF] = get_im_Transforms(im,1);

% Compound transforms: FFT(Gr),FFT(WL),FFT(Bess),FFT(Ch),Ch(Gr),Bess(Gr)
[FFT_Grad,FFT_WL,FFT_Ch,FFT_Bess,Ch_Grad,Bess_Grad] = get_Cmp_Transforms(WL_imTSF,Ch_imTSF,Bess_imTSF,Grad_imTSF);


vec = []; vSizes = [];

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Sigs computed on image...
%[Che,ChF,c4m,edge,fsta,Gab,svds,Har,mhi,rad,tam,zer] = sigs_from_1transform(im);
[foo_vec,foo_sz] = sigs_from_1transform( double(im) );
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Sigs computed on transforms...

% on FFT transform
%[Che_FFT,ChF_FFT,c4m_FFT,edge_FFT,fsta_FFT,Gab_FFT,svds_FFT,Har_FFT,mhi_FFT,rad_FFT,tam_FFT,zer_FFT] = sigs_from_1transform(FFT_imTSF);
[foo_vec,foo_sz] = sigs_from_1transform( double(FFT_imTSF) );
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];

% on Che transform
%[Che_Ch,ChF_Ch,c4m_Ch,edge_Ch,fsta_Ch,Gab_Ch,svds_Ch,Har_Ch,mhi_Ch,rad_Ch,tam_Ch,zer_Ch] = sigs_from_1transform(Ch_imTSF);
[foo_vec,foo_sz] = sigs_from_1transform( double(Ch_imTSF) );
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];

% on Bessel transform
[foo_vec,foo_sz] = sigs_from_1transform( double(Bess_imTSF) );
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];

% on WL transform
%[Che_WL,ChF_WL,c4m_WL,edge_WL,fsta_WL,Gab_WL,svds_WL,Har_WL,mhi_WL,rad_WL,tam_WL,zer_WL] = sigs_from_1transform(WL_imTSF);
[foo_vec,foo_sz] = sigs_from_1transform( double(WL_imTSF) );
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];

% on Gr transform
%[Che_Gr,ChF_Gr,c4m_Gr,edge_Gr,fsta_Gr,Gab_Gr,svds_Gr,Har_Gr,mhi_Gr,rad_Gr,tam_Gr,zer_Gr] = sigs_from_1transform(Grad_imTSF);
[foo_vec,foo_sz] = sigs_from_1transform( double(Grad_imTSF) );
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];

% on Lap transform
%[Che_Lap,ChF_Lap,c4m_Lap,edge_Lap,fsta_Lap,Gab_Lap,svds_Lap,Har_Lap,mhi_Lap,rad_Lap,tam_Lap,zer_Lap] = sigs_from_1transform(Lap_imTSF);
[foo_vec,foo_sz] = sigs_from_1transform( double(Lap_imTSF) );
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Sigs computed on compound transforms...

% on FFT(Gr) c.transform
%[Che_FFT_Gr,ChF_FFT_Gr,c4m_FFT_Gr,edge_FFT_Gr,fsta_FFT_Gr,Gab_FFT_Gr,svds_FFT_Gr,Har_FFT_Gr,mhi_FFT_Gr,rad_FFT_Gr,tam_FFT_Gr,zer_FFT_Gr] = sigs_from_1transform(FFT_Grad_imTSF);
[foo_vec,foo_sz] = sigs_from_1transform(FFT_Grad);
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];

% on FFT(WL) c.transform
[foo_vec,foo_sz] = sigs_from_1transform(FFT_WL);
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];

% on FFT(Ch) c.transform
[foo_vec,foo_sz] = sigs_from_1transform(FFT_Ch);
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];

% on FFT(Bess) c.transform
[foo_vec,foo_sz] = sigs_from_1transform(FFT_Bess);
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];

% on Ch(Gr) c.transform
%[Che_Ch_Gr,ChF_Ch_Gr,c4m_Ch_Gr,edge_Ch_Gr,fsta_Ch_Gr,Gab_Ch_Gr,svds_Ch_Gr,Har_Ch_Gr,mhi_Ch_Gr,rad_Ch_Gr,tam_Ch_Gr,zer_Ch_Gr] = sigs_from_1transform(Ch_Grad);
[foo_vec,foo_sz] = sigs_from_1transform(Ch_Grad);
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];

% on Bess(Gr) c.transform
[foo_vec,foo_sz] = sigs_from_1transform(Bess_Grad);
vec = [vec foo_vec]; vSizes = [vSizes foo_sz];

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Form sig vector
signature_vector = vec;
sigSizes = vSizes;

 % it must be a double type...
if ~isa(signature_vector,'double'), signature_vector = double(signature_vector); end
end % end function


% All 12 sig vectors
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [vec,sigSizes] = sigs_from_1transform(im), vec = [];
%function [Che,ChF,c4m,edge,fsta,Gab,svds,Har,mhi,Rad,Tam,zer] = sigs_from_1transform(im),

if isfloat(im) & ~isa(im,'single'), im = single(im);end
if ~isfloat(im),im = single(im);end
warning off;

% 01: Chebyshev polynomials: 32 elements
try, Che = concat_outputs (ChebyshevStatistics(im));
catch, Che = zeros(1,32);
end
% 02: Ch-F polynomials: 32 elements 
try, ChF = concat_outputs (ChebyshevFourierTransform(im));
catch, ChF = uint16(zeros(1,32));
end
% 03: Comb-4-moments sigs: 48 elements
try, c4m = concat_outputs (CombFirst4Moments(im));
catch, c4m = uint16(zeros(1,48));
end

% 04: Edge sigs
[EdgeArea,MagMean,MagMedian,MagVar,MagHist, DirecMean, DirecMedian,DirecVar,DirecHist,DirecHomo,DiffDirecHist] = EdgeStatistics(double(single(ImageGradient(im))));
edge = ro_EdgeStatistics(EdgeArea,MagMean,MagMedian,MagVar,MagHist, DirecMean, DirecMedian,DirecVar,DirecHist,DirecHomo,DiffDirecHist);

% 05: Feature statistics
[level, BinaryMask] = OtsuGlobalThreshold(double(im));
[Count, Euler, Centroid, AreaMin, AreaMax, AreaMean, AreaMedian, AreaVar, AreaHist, DistMin, DistMax, DistMean, DistMedian, DistVar, DistHist] = ObjectStatistics(logical(BinaryMask));
fsta = ro_FeatureStatistics(Count, Euler, Centroid, AreaMin, AreaMax, AreaMean, AreaMedian, AreaVar, AreaHist, DistMin, DistMax, DistMean, DistMedian, DistVar, DistHist);

% 06: Gabor textures: 28 elements
try, Gab = concat_outputs (GaborTextureFilters2(im));
catch, Gab = zeros(1,28);
end
% 07: SVD sig: 32 elements
try, svds = concat_outputs (SVD_signatures(im));
catch, svds = zeros(1,32);
end
% 08: Haralik series: 28 elements
try, Har = concat_outputs (HaralickTexturesRI(im2uint8(im)));
catch, Har = zeros(1,28);
end
% 09: Multi-scale histogram: 24 elements
try, mhi = concat_outputs (MultiScaleHistograms(im));
catch, mhi = zeros(1,24);
end
% 10: Radon sig: 12 elements
try, Rad = concat_outputs (RadonTransform(double(im)));
catch, Rad = zeros(1,12);
end
% 11: Tamura sig: 6 elements
try, Tam = concat_outputs (TamuraTextures(im));
catch, Tam = zeros(1,6);
end
% 12: Zernike polynomials: 72 elements
try, Zer = concat_outputs (mb_zernike(im));
catch, Zer = zeros(1,72);
end

vec = [Che ChF c4m edge fsta Gab svds Har mhi Rad Tam Zer];
sigSizes = [length(Che) length(ChF) length(c4m) length(edge) length(fsta) length(Gab) length(svds) length(Har) length(mhi) length(Rad) length(Tam) length(Zer)];
end % eofunc


% Compound transforms: FFT(Gr),FFT(WL),FFT(Bess),FFT(Ch),Ch(Gr),Bess(Gr)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [FFT_Grad,FFT_WL,FFT_Ch,FFT_Bess,Ch_Grad,Bess_Grad] = get_Cmp_Transforms(WL_imTSF,Ch_imTSF,Bess_imTSF,Grad_imTSF),
%function [FFT_Grad,FFT_WL,FFT_Ch,FFT_Bess,Ch_Grad,Bess_Grad] = get_Cmp_Transforms(WL_imTSF,Ch_imTSF,Bess_imTSF,Grad_imTSF),

if ~isa(WL_imTSF,'double'),WL_imTSF = double(WL_imTSF);end
if ~isa(Ch_imTSF,'double'),Ch_imTSF = double(Ch_imTSF);end
if ~isa(Bess_imTSF,'double'),Bess_imTSF = double(Bess_imTSF);end
if ~isa(Grad_imTSF,'double'),Grad_imTSF = double(Grad_imTSF);end

FFT_Grad  = abs( fftshift(fft2( Grad_imTSF )) ); FFT_Grad = log( FFT_Grad ); FFT_Grad = norm0max( FFT_Grad );
FFT_WL    = abs( fftshift(fft2( WL_imTSF )) );   FFT_WL = log( FFT_WL );     FFT_WL = norm0max( FFT_WL );
FFT_Ch    = abs( fftshift(fft2( Ch_imTSF )) );   FFT_Ch = log( FFT_Ch );     FFT_Ch = norm0max( FFT_Ch );
FFT_Bess  = abs( fftshift(fft2( Bess_imTSF )) ); FFT_Bess = log( FFT_Bess ); FFT_Bess = norm0max( FFT_Bess );

Ch_Grad   = ChebyshevTransform(Grad_imTSF);
Bess_Grad = BesselTransform( Grad_imTSF );

%vec = [FFT_Grad FFT_WL FFT_Ch FFT_Bess Ch_Grad Bess_Grad];
%sigSizes = [length(FFT_Grad) length(FFT_WL) length(FFT_Ch) length(FFT_Bess) length(Ch_Grad) length(Bess_Grad)];
end % eofunc

% Computes Transforms on image
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [FFT_imTSF,Ch_imTSF,Bess_imTSF,WL_imTSF,Grad_imTSF,Lap_imTSF] = get_im_Transforms(im,level),
%if isfloat(im) & ~isa(im,'single'), im = single(im);end
%if ~isfloat(im),im = single(im);end
% FFT(im)
%FFT_imTSF = FrequencySpace2Pixels(FourierTransform(im));
FFT_imTSF = abs( fftshift(fft2( im )) ); FFT_imTSF = log( FFT_imTSF ); FFT_imTSF = norm0max( FFT_imTSF );
% Cheb(im)
Ch_imTSF = ChebyshevTransform(im); %Ch_imTSF = double(Ch_imTSF);
% Bessel(img)
Bess_imTSF = BesselTransform( im );
% WL(im)
WL_imTSF = WaveletSelector(WaveletTransform(im)); %WL_imTSF = double(WL_imTSF);

% Gradient(im) and Laplacian(im)
Grad_imTSF = []; Lap_imTSF = [];
if level==1, [Grad_imTSF,Lap_imTSF] = get_diff_transforms(im);end % level
end % eofunc


% Computes Gradient and Laplacian of the image
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [grad,Lap,gradgr] = get_diff_transforms(Im),
% determine type of image array
if isfloat(Im),
 if isa(Im,'double'),imtype = 'double';
 else, imtype = 'single'; Im = double(Im);
 end
else,
 if isa(Im,'uint8'),imtype = 'uint8';end
 if isa(Im,'uint16'),imtype = 'uint16';end
 Im = double(Im);
end % if isfloat(Im)
% Gradient
[Fx,Fy] = gradient(Im); grad = norm01(sqrt(Fx.^2 + Fy.^2));
% Laplacian
[Fxx,Fxy] = gradient(Fx); [Fyx,Fyy] = gradient(Fy);
Lap = norm01(Fxx.^2 + Fyy.^2);
return;
switch imtype,
 case 'uint8', Lap = uint8(255.*Lap); grad = uint8(255.*grad);
 case 'uint16', Lap = uint16(65535.*Lap); grad = uint16(65535.*grad);
 case 'single', Lap = single(Lap); grad = single(grad);
end % switch
end % eofunc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function y = norm0max( x ), x0 = min(x(:)); x1 = max(x(:));
y = x1.* (x-x0)./(x1-x0+eps);
end % eofunc
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function y = norm01( x ), x0 = min(x(:)); x1 = max(x(:));
y = (x-x0)./(x1-x0+eps);
end % eofunc


%
% MATLAB implementation of required typecaster modules
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function pix = FrequencySpace2Pixels(fs)
pix = single(fs(:,:,1)); % the single is cause OMEIS doesn't support doubles
end % end function
 
function wav_result = WaveletSelector(wav1, wav2)
wav_result = wav1; % the single is cause OMEIS doesn't support doubles
end % end function

%
% MATLAB implementation of certain module's vector decoder (from execution instructions)
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function out = vd_Comb4Moments(in)
out = [in(46) in(47) in(48) in(37) in(38) in(39) in(43) in(44) in(45) in(40) in(41) in(42) ...
   in(34) in(35) in(36) in(25) in(26) in(27) in(31) in(32) in(33) in(28) in(29) in(30) ...
   in(10) in(11) in(12) in(1) in(2) in(3) in(7) in(8) in(9) in(4) in(5) in(6) ...
   in(22) in(23) in(24) in(13) in(14) in(15) in(19) in(20) in(21) in(16) in(17) in(18)];
end % end function
   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function out = ro_EdgeStatistics(EdgeArea, MagMean, MagMedian, MagVar, MagHist, DirecMean, DirecMedian, DirecVar, DirecHist, DirecHomogeneity, DiffDirecHist)
out = [double(EdgeArea) double(DiffDirecHist) double(DirecHist) double(DirecHomogeneity)...
	   double(DirecMean) double(DirecMedian) double(DirecVar) double(MagHist) double(MagMean) double(MagMedian) double(MagVar)];
end % end function

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function out = ro_FeatureStatistics(Count, Euler, Centroid, AreaMin, AreaMax, AreaMean, AreaMedian, AreaVar, AreaHist, DistMin, DistMax, DistMean, DistMedian, DistVar, DistHist)
out = [double(AreaHist) double(AreaMax) double(AreaMean) double(AreaMedian) ...
	   double(AreaMin) double(AreaVar) double(Centroid) double(Count) double(DistHist) ...
	   double(DistMax) double(DistMean) double(DistMedian) double(DistMin) double(DistVar) double(Euler)];
end % end function
	   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function out = vd_HaralickTexturesRI(in)
out = [in(1) in(15) in(2) in(16) in(3) in(17) in(10) in(24) in(11) in(25) in(9) in(23) ...
   in(12) in(26) in(5) in(19) in(14) in(28) in(13) in(27) in(6) in(20) in(8) in(22) ...
   in(7) in(21) in(4) in(18)];
end % end function

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function out = vd_RadonTextures(in)
out = [in(1) in(2) in(3) in(10) in(11) in(12) in(4) in(5) in(6) in(7) in(8) in(9)];
end % end function

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function out = vd_TamuraTextures(in)
out = [in(2) in(3) in(4) in(6) in(5) in(1)];
end % end function

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
