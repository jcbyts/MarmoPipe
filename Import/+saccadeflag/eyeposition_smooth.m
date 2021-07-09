function [smo] = eyeposition_smooth( eyePos , varargin)
% SMOOTH EYE POSITION
% [smo] = eyeposition_smooth( eye )
%
%
% Input:
%  eyePos [T x 4] where T is number of samples
%       column 1: time (seconds)
%       column 2: x position (dva)
%       column 3: y position (dva)
%       column 4: pupil size
%
% Output:
%  eyeSmo [Tnew x 4] similar data struct, but up-sampled to 1kHz
%       column 1: time (seconds)
%       column 2: x position (dva)
%       column 3: y position (dva)
%       column 4: pupil size
%       column 5: x velocity
%       column 6: y velocity
%       column 7: eye speed
%
% Optional arguments:
%
% 'Smoothing' (10) Gaussian smoothing (in ms)
%
%*******************************************

ip = inputParser();
ip.KeepUnmatched = true; % ignore unmatched input args
ip.addParameter('Smoothing', 2);  %this was 10 or 5 ms in Sunwoo's
ip.parse(varargin{:});

Smoothing = ip.Results.Smoothing;

%****** use linear interpolate to up-sampled to 1 kHz
tt = eyePos(:,1);
[~,xx2] = interp_data(tt,eyePos(:,2));
[~,yy2] = interp_data(tt,eyePos(:,3));
[tt2,pp2] = interp_data(tt,eyePos(:,4));
%****** use Gauss smoothing to smooth traces
[~,xx3] = smooth_data(tt2,xx2,Smoothing);
[~,yy3] = smooth_data(tt2,yy2,Smoothing);
[tt3,pp3] = smooth_data(tt2,pp2,Smoothing);
%******* differentiate for velocity
[vx,vy] = differentiate(tt3,xx3,yy3);
av = hypot(vx, vy); % fast and numerically stable hypotenuse
%******** return smoothed and velocity data
smo = [tt3(:), xx3(:), yy3(:), pp3(:), vx(:), vy(:), av(:)];

return;

function [tt2,xx2] = interp_data(tt,xx)
%*********
if (max(size(tt)) < 10)
    tt2 = NaN;
    xx2 = NaN;
    return
end
%************
tt2 = tt(1):0.001:(max(tt)-0.001);  % resampled data up to max time (1000 kHz)
npoints = length(tt2);
xx2 = nan(1,npoints);

for zk = 1:npoints
    targtime = tt2(zk);
    %********** find where it is in between
    for kk = 1:(size(tt,1)-1)  % go to actual, search for where new time is located between two points
        if (targtime >= tt(kk)) && (targtime < tt(kk+1))
            break;
        end
    end
    %** now at kk you are in between the two points
    XA = xx(kk);   % lower value before newtime
    XB = xx(kk+1); % upper value after newtime
    t1 = targtime - tt(kk);
    t2 = tt(kk+1) - targtime;
    %*********
    val = (((t2/(t1+t2))*XA) + (t1/(t1+t2))*XB);
    xx2(zk) = val;
    %**********
end


function [tt3,xx3] = smooth_data(tt2,xx2,GSig)
pts = GSig; % gaussian size
dpts = (2*GSig)+1;
posfilt = gengauss(pts,dpts); % create a gaussian filter (20 samples)
posfilt = posfilt/sum(posfilt);
xx3 = conv(posfilt,xx2);
xx3 = xx3(pts:end-(pts+1)); % Window back
xx3(1:pts) = xx3(pts+1);  % flat position on leading edge (ploblem of filtering)
xx3((end-pts+1):end) = xx2(end-pts);  % again flat
tt3 = tt2;


function x = gengauss(w,p,c)
%GENGAUSS Generate gaussian peak of unity height.
%   X = GENGAUSS(W)  or   X = GENGAUSS(W,P,C)
%    where:
%    W    is the distance from center to 2% of full height
%    P    (optional) is the number of points x will contain
%    C    (optional) is the point number for the peak center
%    The peak will be centered in the vector unless P and C are specified.
%   by Richard Kramer.
%   Copyright (c) 1988-1993 by The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 1993/09/01 19:37:08 $

if nargin == 2, c = 0; end
if nargin == 1, c = 0; p = 0; end

if  p < 1
    points = round(3 * w);
else
    points = round(p);
end

if c < 1
    center = round(points/2);
else
    center = round(c);
end

w = w/2;
if c < 0, center = center - c; end

x = nan(1,points);
for n = 1:points
    x(n) = exp(-((n-center)/w)^2);
end


%******** computes the derivative
function dvec = differ(mat)
dvec = zeros(size(mat));
for k = 1:length(dvec)
    if (k == 1)
        dvec(k) = 2*(mat(k+1)-mat(k));
    else
        if (k == length(mat))
            dvec(k) = 2*(mat(k)-mat(k-1));
        else
            dvec(k) = (mat(k+1)-mat(k-1));
        end
    end
end
dvec = dvec * 1000;  % go from milliseconds to seconds


%compute derivative numerically of a trace
function [vx,vy] = differentiate(tt3,xx3,yy3)
vx = (differ(xx3)/2)';
vy = (differ(yy3)/2)';

