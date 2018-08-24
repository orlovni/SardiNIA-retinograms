%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function tvec = tortuosity_15vec ( xraw,yraw ), tvec = [];

% Coordinates must be column vectors
xraw = xraw(:);
yraw = yraw(:);

% desired segment length, sampling parameter
desiredSL = 9;

%[xs,ys] = xy_discretize(xraw,yraw,desiredSL);
xs = xraw; ys = yraw;
[DF,t2,t3,t4,t5,t6,t7,t8] = tortuosity_lib( xs,ys,desiredSL );

tvec = DF(:);
tvec = [tvec; t2(:)];
tvec = [tvec; t3(:)];
tvec = [tvec; t4(:)];
tvec = [tvec; t5(:)];
tvec = [tvec; t6(:)];
tvec = [tvec; t7(:)];
tvec = [tvec; t8(:)];

end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function varargout = tortuosity_lib( xr,yr,dSL ),
[xi,yi] = xy_spline( xr,yr,dSL );
xi = xi(:);
yi = yi(:);

% Normalized arc length for sampled and interpolated vessel midline
ALs = arcLen( xr,yr );
CLs = cord_len( xr,yr );
DF1 = ALs ./ CLs;

ALi = arcLen( xi,yi );
CLi = cord_len( xi,yi );
DF2 = ALi ./ CLi;

varargout{1} = [DF1 DF2];

% Grad and GradGrad for sampled and interpolated vessel midline
x1s = gradient(xr); y1s = gradient(yr);
x2s = gradient(x1s);  y2s = gradient(y1s);
x1i = gradient(xi);   y1i = gradient(yi);
x2i = gradient(x1i);  y2i = gradient(y1i);

% Total curvature for sampled and interpolated vessel midline
tau2_sampled = totalCurvature( x1s,y1s,x2s,y2s );
tau2_interpd = totalCurvature( x1i,y1i,x2i,y2i );
varargout{2} = [tau2_sampled tau2_interpd];

% Total squared curvature for sampled and interpolated vessel midline
tau3_sampled = totalSquaredCurvature( x1s,y1s,x2s,y2s );
tau3_interpd = totalSquaredCurvature( x1i,y1i,x2i,y2i );
varargout{3} = [tau3_sampled tau3_interpd];

% Total curvature / arc length
tau4_sampled = tau2_sampled ./ ALs;
tau4_interpd = tau2_interpd ./ ALi;
varargout{4} = [tau4_sampled tau4_interpd];

% Total squared curvature / arc length
tau5_sampled = tau3_sampled ./ ALs;
tau5_interpd = tau3_interpd ./ ALi;
varargout{5} = [tau5_sampled tau5_interpd];


% Total curvature / chord length
tau6_sampled = tau2_sampled ./ CLs;
tau6_interpd = tau2_interpd ./ CLi;
varargout{6} = [tau6_sampled tau6_interpd];

% Total squared curvature / chord length
tau7_sampled = tau3_sampled ./ CLs;
tau7_interpd = tau3_interpd ./ CLi;
varargout{7} = [tau7_sampled tau7_interpd];

% Sampled arc length / interpolated arc length
tau10 = ALi / ALs;
varargout{8} = tau10;
end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function TC = totalCurvature( x1,y1,x2,y2 )
%
% Total curvature for vessel midline
% Bhuiyan et al (2010)
% Hart et al (1997,1999)
%
foo = x1.*y2 - x2.*y1;
arcLen = ( x1.^2 + y1.^2 ).^(1.5);
TC = sum( foo./arcLen );
end % eofunc



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function AL = arcLen( x,y )
%
% Normalized arc length for vessel midline.
% Benitez-Aguire et al (2011,2012)
% Zepeda-Romero et al (2011)
% Mahal et al (2009)
%
AL = arc_len( x,y );
end % eofunc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function arcLen = discrete_archLen( x,y )
x1 = gradient(x); y1 = gradient(y);
arcLen = sqrt( x1.^2 + y1.^2 );
end % eofunc




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function TSC = totalSquaredCurvature( x1,y1,x2,y2 )
%
% Total squared curvature for vessel midline
% Hart et al (1997,1999)
%
foo = x1.*y2 - x2.*y1;
arcLen = ( x1.^2 + y1.^2 ).^(1.5);
TSC = sum( foo.^2./arcLen.^2 );
end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function AL = arc_len( x,y )
%
% Arc length for vessel midline.
%
x1 = gradient(x); y1 = gradient(y);
AL = sqrt( x1.^2 + y1.^2 );
AL = sum( AL );
end % eofunc

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function CL = cord_len( x,y )
%
% Cord length for vessel midline.
%
v1 = [x(1) - x(end); y(1) - y(end)];
CL = norm( v1 );
end % eofunc




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [x2,y2,xAnchors,yAnchors] = xy_spline( xr,yr,dSL )
%
% To smooth the curve (x,y) with spline
%
n = length(xr);
if dSL > 2*n,
 x2 = xr(:)'; y2 = yr(:)';
 xAnchors = xr(:)'; yAnchors = yr(:)';
 return;
end

abscissa = 1:n;
indAnchors = 1:round(dSL):n;
if (indAnchors(end) < n) & (n - indAnchors(end) > 3),indAnchors = [indAnchors n];end
xAnchors = xr(indAnchors);
yAnchors = yr(indAnchors);
x2 = interp1( indAnchors,xAnchors,abscissa,'spline' );
y2 = interp1( indAnchors,yAnchors,abscissa,'spline' );
end % eofunc


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [x,y,x2,y2,xAnchors,yAnchors] = xy_discretize( xr,yr,dSL )
%
% To interpolate the curve (x,y) in a piece-wise linear fashion
%
n = length(xr);
x = []; y = [];
x2 = []; y2 = [];

% One way I may do, is to go by intersects of the curve and the circle of
% radius dSL. E.g. the cost function is CF = (xx-x0).^2 + (yy-y0).^2 - dSL.^2
% However, if the curve is too windy I get multiple intersections of the
% circle and the curve -- no good.
% Another way is more simple. The curve parameter (new abscissa) is the
% indices of (xr,yr); it is absolutely linear.

abscissa = 1:n;
indAnchors = 1:round(dSL):n;
xAnchors = xr(indAnchors);
yAnchors = yr(indAnchors);
x  = interp1( indAnchors,xAnchors,abscissa(1:indAnchors(end)),'linear' );
y  = interp1( indAnchors,yAnchors,abscissa(1:indAnchors(end)),'linear' );
x2 = interp1( indAnchors,xAnchors,abscissa(1:indAnchors(end)),'PCHIP' );
y2 = interp1( indAnchors,yAnchors,abscissa(1:indAnchors(end)),'PCHIP' );

% The curve (x2,y2) is not meant indices; it is rather for interim
% computation of tortuosity.
%x2 = round(x2);y2 = round(y2);
end % eofunc


