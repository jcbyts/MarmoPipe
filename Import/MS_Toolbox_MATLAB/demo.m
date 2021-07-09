% % #-------------------------------------------------
% % # demonstration for microsaccade detection and
% % # surrogate analysis
% % # Ralf Engbert (2013)
% % #-------------------------------------------------
clear all

% set parameter
SAMPLING = 500;
MINDUR = 6*SAMPLING/1000;
VFAC = 5;
% % read raw data
d = load(sprintf('data/f0%d.00%d.dat',1 , 4));
% % select epoch, transform to matrix
idx = 3001:4500;
xl = d(idx,2:3);
xr = d(idx,4:5);

% % #-------------------------------------------------
% % # detect microsaccades 
% % #-------------------------------------------------
% % detect microsaccades
[msl, radiusl] = microsacc(xl,VFAC,MINDUR,SAMPLING);
msr = microsacc(xr,VFAC,MINDUR,SAMPLING);
sac = binsacc(msl,msr);
N = length(sac(:,1));
xls = smoothdata(xl);

% figure(1),clf
% plot(xl(:,1),'b.'),hold on
% plot(xl(:,2),'r.'),hold on
% xlabel('Time (ms)')
% ylabel('Positions')
% legend('Horizontal','Vertical')
% 
% x_smooth = xls(:,1);
% y_smooth = xls(:,2);
% 
% plot(x_smooth,'b-')
% plot(y_smooth,'r-')
% 
% pause

% % plot trajectory
h=figure(1);set(h,'name','Microsaccade detection')
subplot(1,2,1);
plot(xls(:,1),xls(:,2),'k'); hold on;
xlabel('xl');ylabel('yl'); title('Position');
for s = 1:N 
    j = sac(s,1):sac(s,2); 
    plot(xls(j,1),xls(j,2),'r');
end
plot(xls(sac(:,2),1),xls(sac(:,2),2),'or')
hold off; axis square; 

% % plot velocity data
subplot(1,2,2)
vls = vecvel(xl,SAMPLING);
plot(vls(:,1),vls(:,2),'k');hold on;
xlabel('v(x)');ylabel('v(y)'); title('Velocity');
for s = 1:N 
    j = sac(s,1):sac(s,2); 
    plot(vls(j,1),vls(j,2),'r')
end
phi = linspace(0,2*pi,300);
cx = radiusl(1)*cos(phi);
cy = radiusl(2)*sin(phi);
plot(cx,cy,':k');
hold off; axis square; 


% % #-------------------------------------------------
% % # generate surrogates and polt acf
% % #-------------------------------------------------
% % generate surrogates
xlsur = surrogate(xls,SAMPLING);

% % plot original trajectory, surrogate trajectory
h=figure(2);set(h,'name','Surrogates')
subplot(2,2,1)
plot(xls(:,1),xls(:,2),'k')
hold off; axis square; 
xlabel('x'); ylabel('y'); title('Original data');
subplot(2,2,2)
plot(xlsur(:,1),xlsur(:,2),'r')
hold off; axis square; 
xlabel('x'); ylabel('y'); title('AAFT surrogate');

% % plot acf of original trajectory and surrogate trajectory
subplot(2,2,3)
vls = vecvel(xls,SAMPLING);
[a1,lag] = acorr(vls(:,1));
plot(a1,'k'); hold on;
for s = 1:10
    xlsur = surrogate(xls,SAMPLING);
    vlsur = vecvel(xlsur,SAMPLING);
    a2 = acorr(vlsur(:,1));
    plot(a2,'r');
end
plot(lag(1):lag(end):lag(end),[ 0  0],'k'); 
plot(lag(1):lag(end):lag(end),[.05 .05],'--b'); 
plot(lag(1):lag(end):lag(end),[-.05 -.05],'--b'); 
hold off; axis square;
xlabel('Lag'); ylabel('ACF'); title('Autocorrelation function');
 

% % #-------------------------------------------------
% % # surrogate analysis for microsaccade detection
% % #-------------------------------------------------
vfac = 3:.5:8;  % range of detection thresholds 
vp = 1;         % subjects and trials
ntrials = 5;    % number of trials
mstab = zeros(length(vfac),3);
for v = 1:vp
    for n = 1:ntrials
        dd = load(sprintf('data/f0%d.00%d.dat', v, n));
        xxl = dd(:,2:3);
        xxr = dd(:,4:5);
        dur = length(d(:,1))/SAMPLING;        
        for i = 1:length(vfac)
            % % detect microsaccades
            msl = microsacc(xxl,vfac(i),MINDUR,SAMPLING);
            msr = microsacc(xxr,vfac(i),MINDUR,SAMPLING);
            sac = binsacc(msl, msr);
            if isempty(sac); N = 0;
            else             N = length(sac(:,1))/dur; 
            end
            % % computation with surrogate data
            xlsur = surrogate(xxl,SAMPLING);
            xrsur = surrogate(xxr,SAMPLING);
            msl = microsacc(xlsur,vfac(i),MINDUR,SAMPLING);
            msr = microsacc(xrsur,vfac(i),MINDUR,SAMPLING);
            sac = binsacc(msl, msr);
            if isempty(sac); Nsur = 0;
            else             Nsur = length(sac(:,1))/dur;
            end
            mstab(i,1) = vfac(i);
            mstab(i,2:3) = mstab(i,2:3) + [N Nsur];
        end
    end
end
mstab(:,2:3) = mstab(:,2:3)./(vp*ntrials);

% % plot surrogate analysis
h=figure(3);set(h,'name','Surrogate analysis')
plot(mstab(:,1),mstab(:,2),'o-k'); hold on;
plot(mstab(:,1),mstab(:,3),'o-r');
plot(mstab(:,1),mstab(:,2)-mstab(:,3),'o-b');
legend('original data', 'surrogates', 'difference');
xlabel('\lambda'); ylabel('Rate [1/s]'); hold off;


% % #-------------------------------------------------
% % # random walk analysis
% % #-------------------------------------------------
Np = length(xl(:,1));
rvl = lagdist(xls);
vl = diff(xls);
vls = [vl(randperm(Np),1), vl(randperm(Np),2)];
rxls = [cumsum(vls(:,1)),cumsum(vls(:,2))];
rvls = lagdist(rxls);

% % add plot for simulated data (SAW model)
ds = load('data/output.dat');
sim = ds(:,2:3);
rsim = lagdist(sim);
rsim(2,:) = rsim(2,:)./5000;

% % plot random walk analysis
h=figure(4);set(h,'name','Random walk analysis')
loglog(rvl(1,:),rvl(2,:),'k'); hold on;
loglog(rvls(1,:),rvls(2,:),'r');
loglog(rsim(1,:),rsim(2,:),'b')
hold off;  
xlabel('Lag');ylabel('D(lag)'); title('Random walk analysis');


% % #-------------------------------------------------
% % # box-count analysis
% % #-------------------------------------------------
dx = 0.01;          % edge length dx
dt = 100;           % time window
msl = microsacc(xls,VFAC,MINDUR,SAMPLING);

% % calculate box-count of drift periode
xx = xls(msl(1,1)-dt:msl(1,1)-1,:);
boxes = boxcount(xx,dx);

% % plot box counting of drift periode
h=figure(5);set(h,'name','Box-count')
plot(xx(:,1),xx(:,2));hold on;
x_min = min(xx(:,1));
y_min = min(xx(:,2));
M = length(xx);
for l=1:M
    i = floor( (xx(l,1)-x_min)/dx ) + 1;
    j = floor( (xx(l,2)-y_min)/dx ) + 1;
    plot([x_min+i*dx-dx,x_min+i*dx],[y_min+j*dx-dx,y_min+j*dx-dx],'g')
    plot([x_min+i*dx-dx,x_min+i*dx],[y_min+j*dx,y_min+j*dx],'g')
    plot([x_min+i*dx-dx,x_min+i*dx-dx],[y_min+j*dx-dx,y_min+j*dx],'g')
    plot([x_min+i*dx,x_min+i*dx],[y_min+j*dx-dx,y_min+j*dx],'g')
end
text(max(xx(:,1))-dx*4,min(xx(:,2))+dx,sprintf('d=%i in %i ms',boxes, dt));
xlabel('x [deg]');ylabel('y [deg]'); title('Box-count');
hold off;
