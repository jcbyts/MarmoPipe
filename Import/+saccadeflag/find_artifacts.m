function artifacts = find_artifacts(timestamps, xpos, ypos,varargin)
% Find saccades in eye position signals...
%
% Available arguments include:
%
%   order     - low-pass digital differentiating filter order (default: 32)
%   Wn        - low-pass filter corner freq, and transition band as a
%               percentage of the Nyquist frequency (default: [0.1,0.16])
%   accthresh - acceleration threshold (default: 2e4 deg./s^2)
%   velthresh - velocity threshold (default: 10 deg./s)
%   velpeak   - minimum peak velocity (default: 10 deg./s)
%   isi       - minimum inter-saccade interval (default: 0.050s)
%   debug     - show debugging output (default: false)
%
% See also: lpfirdd.

% 2023-06-20 - Jake

args = varargin;
p = inputParser;

% p.addParameter('order',5,@(x) validateattributes(x,{'numeric'},{'scalar','even'})); % low-pass filter/differentiator order
% p.addParameter('Wn',[0.1,0.16],@(x) validateattributes(x,{'numeric'},{'vector'})); % filter transition band (percentage of Nyquist frequency)
p.addParameter('order',5);
p.addParameter('velthresh',10,@(x) validateattributes(x,{'numeric'},{'scalar','positive'})); % velocity threshold (deg./s)
p.addParameter('buffer', 0)
p.addParameter('isi',0.050,@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
p.addParameter('dt',0.075,@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
p.addParameter('debug',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

p.parse(args{:});

args = p.Results;

% low-pass FIR digital differentiator coefficients
% N = round(args.order/2);
% copt = lpfirdd(N, args.Wn(1), args.Wn(2), 1 ,0)';
% coeffs = [fliplr(copt), 0, -copt];

timestamps = timestamps(:);
xpos = xpos(:);
ypos = ypos(:);

% get gaze position data...
t = timestamps;
pos = struct('x',xpos,'y',ypos);

fs = round(1/nanmedian(diff(timestamps))); % sampling freq. (samples/s)
% horiz. (x) and vert. (y) velocities
dxdt = @(x) imgaussfilt(filter([1; -1], 1, x), args.order);
velx = dxdt(xpos)*fs;
vely = dxdt(ypos)*fs;

% vel = structfun(@(x) fs*locfilt(coeffs,1,x),pos,'UniformOutput',false);


% vel = structfun(@(x) fs*dxdt(x),pos,'UniformOutput',false);
% scalar eye speed
speed = hypot(velx,vely);
speed = speed + 1./speed;

% % estimate baseline (e.g., pursuit) speed using a moving average...
% a = 1;
% b = ones(1,2*N+1)./(2*N+1);
% baseline = locfilt(b,a,speed);

% find baseline velocities that exceed threshold
startidx = findZeroCrossings(fix(speed - args.velthresh),1);
stopidx = findZeroCrossings(fix(speed - args.velthresh),-1);

a = 1;

if stopidx(end) < startidx(end)
    stopidx = [stopidx; numel(speed)];
end

if startidx(1) > stopidx(1)
    startidx = [1; startidx];
end

while numel(startidx) > numel(stopidx)
    
    if (stopidx(a) - startidx(a+1)) >= 0
        startidx(a + 1) = [];
    end
    a = a + 1;
end

% This has not been debugged
a = 1;
while numel(startidx) < numel(stopidx)
    
    if (stopidx(a) - startidx(a)) <= 0
        stopidx(a) = [];
    end
    a = a + 1;
end


%%

if numel(startidx) > 1

    if startidx(1) > stopidx(1)
        startidx = [1; startidx];
    end
    
    if stopidx(end) < startidx(end)
        stopidx = [stopidx; numel(baseline)];
    end

    assert(numel(stopidx)==numel(startidx), 'mismatch in size')

    n = ceil(args.isi*fs); % samples
    
    % test for minimum inter-block interval (args.isi) violations?
    short = find(startidx(2:end) - stopidx(1:end-1) < n);
    startidx(short+1) = [];
    stopidx(short) = [];
    

end

if args.buffer > 0
    fprintf('adding buffer %d\n', args.buffer)
    startidx = max(1, startidx - args.buffer);
    stopidx = min(numel(t), stopidx + args.buffer);
end

if args.debug
    figure(99); clf
    subplot(2,1,1)
    plot(t, xpos)
    hold on
    plot(t, ypos)
    for i = 1:numel(startidx)
        fill(t([startidx(i), startidx(i), stopidx(i), stopidx(i)]), [ylim fliplr(ylim)]', 'k', 'FaceAlpha', .25, 'EdgeColor','none');
    end
    subplot(2,1,2)
    plot(t, speed); hold on
    for i = 1:numel(startidx)
        fill(t([startidx(i), startidx(i), stopidx(i), stopidx(i)]), [ylim fliplr(ylim)]', 'k', 'FaceAlpha', .25, 'EdgeColor','none');
    end
end



tstart = t(startidx);
tend = t(stopidx);

artifacts = struct();
artifacts.tstart = tstart(:);
artifacts.tend   = tend(:);
artifacts.duration = artifacts.tend- artifacts.tstart;
artifacts.startIndex = startidx;
artifacts.endIndex   = stopidx;


function y = locfilt(b,a,x)
    y = filter(b,a,x);
    y(1:2*N) = NaN;
    y = circshift(y,-N);
end



%---------------------------------------------------------
% Low-pass FIR digital differentiator (LPFIRDD) design   -
% via constrained quadratic programming (QP)             -
% [Copt,c]=lpfirdd(N,alpha,beta,r,idraw)                 -
% By Dr Yangquan Chen		019-07-1999                   -
% Email=<yqchen@ieee.org>; URL=http://www.crosswinds.net/~yqchen/
% --------------------------------------------------------
% LPFIRDD: only 1st order derivative estimate
% total taps=2N. c(i)=-c(i+N+1); c(N+1)=0 (central point)
%
%          -N      -N+1            -1                    N
% FIR=c(1)z   +c(2)z    +...+ c(N)z  + 0 + ... + c(2N+1)z
%
%       N
%     ------
%     \                    j   -j
%      >       Copt(j) * (z - z  )
%     /
%     ------
%      j=1
% N: Taps  (N=2, similar to sgfilter(2,2,1,1)
% alpha ~ beta: transit band of frequency
%				    (in percentage of Nyquest freq)
% r: the polynomial order. r<=N Normally, set it to 1.
%---------------------------------------------------------------
    function [Copt,bd]=lpfirdd(N,alpha,beta,r,idraw)
        % testing parameters
        % alpha=1./pi;beta=1.5/pi;N=10;r=1;idraw=1;
        if (alpha>beta)
            disp('Error in alpha (alpha<=beta)');return;
        end
        if ((beta>1) || (beta <0))
            disp('Error in Beta! (beta in [0,1]');return;
        end
        if ((alpha>1) || (alpha <0))
            disp('Error in Alpha! (Alpha in [0,1]');return;
        end
        % default r=1
        if (r<1); r=1; end
        
        % matrix W
        W=zeros(r,N);
        for ix=1:N
            for jx=1:r
                W(jx,ix)=ix^(2*jx-1);
            end
        end
        
        %matrix L
        L=zeros(N,1);
        if (beta>alpha)
            for ix=1:N
                L(ix)=(alpha*sin(ix*beta*pi)-beta*sin(ix*alpha*pi))/ix/ix/(beta-alpha);
            end
        elseif (beta==alpha)
            for ix=1:N
                L(ix)=(ix*alpha*pi*cos(ix*alpha*pi)-sin(ix*alpha*pi))/ix/ix;
            end
        end
        % matrix e
        ex=zeros(r,1);ex(1)=1;
        % optimal solution
        % Copt=W'*inv(W*W')*(ex + 2.*W*L/pi)-2.*L/pi;
        Copt=W'*pinv(W*W')*(ex + 2.*W*L/pi)-2.*L/pi;
        Copt=Copt/2;
        % fr plots
        if (idraw==1)
            bd=[-fliplr(Copt'),0,Copt']';
            %ad=1;sys_sg=tf(bd',ad,1./Fs);bode(sys_sg)
            Fs=12790;nL=N;nR=N;npts=1000;%w=logspace(0,4,npts);
            w=((1:npts)-1)*pi/npts;
            j=sqrt(-1);ejw=zeros(nL+nR+1,npts);
            for ix=(-nL:nR)
                ejw(ix+nL+1,:)=exp(j*ix*w);
            end
            freq=bd'*ejw;
            figure;subplot(2,1,1)
            plot(w/pi*Fs/2,(abs(freq)));grid on;
            hold on; ax=axis;ax(2)=Fs/2;axis(ax);
            xlabel('freq. (Hz)');ylabel('amplitude (dB)');
            subplot(2,1,2)
            plot(w/pi*Fs/2,180*(angle(freq))/pi );grid on;
            hold on; ax=axis;ax(2)=Fs/2;axis(ax);
            xlabel('freq. (Hz)');ylabel('phase anlge (deg.)');
            
            figure;subplot(2,1,1);Fs=12600; % Hz for U8
            semilogx(w*Fs/pi/2,20*log10(abs(freq)));grid on;
            hold on; ax=axis;ax(2)=Fs/2;axis(ax);
            semilogx([Fs/2,Fs/2],[ax(3),ax(4)],'o-r');grid on;
            xlabel('freq. (Hz)');ylabel('amplitude (dB)');
            subplot(2,1,2)
            semilogx(w*Fs/pi/2,180*(angle(freq))/pi );grid on;
            hold on; ax=axis;ax(2)=Fs/2;axis(ax);
            semilogx([Fs/2,Fs/2],[ax(3),ax(4)],'o-r');grid on;
            xlabel('freq. (Hz)');ylabel('phase anlge (deg.)');
        end
    end

    function indices = findZeroCrossings(data, mode)
        %FINDZEROCROSSINGS Find zero crossing points.
        %   I = FINDZEROCROSSINGS(DATA,MODE) returns the indicies into the supplied
        %   DATA vector, corresponding to the zero crossings.
        %
        %   MODE specifies the type of crossing required:
        %     MODE < 0 - results in indicies for the -ve going zero crossings,
        %     MODE = 0 - results in indicies for ALL zero crossings (default), and
        %     MODE > 0 - results in indicies for the +ve going zero crossings.
        
        % $Id: findZeroCrossings.m,v 1.1 2008-07-21 23:31:50 shaunc Exp $
        
        if nargin < 2
            mode = 0;
        end
        
        [indices,~,p0] = find(data); % ignore zeros in the data vector
        
        switch sign(mode)
            case -1
                % find -ve going crossings
                iit = find(diff(sign(p0))==-2);
            case 0
                % find all zero crossings
                iit = find(abs(diff(sign(p0)))==2);
            case 1
                % find +ve going crossings
                iit = find(diff(sign(p0))==2);
        end
        
        indices = round((indices(iit)+indices(iit+1))/2);
    end

end




