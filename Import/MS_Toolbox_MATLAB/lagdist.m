function rv = lagdist(x)
%---------------------------------------------------------------------
% Compute lagged squared-distance estimator
% Random walk analysis
% 2012, Ralf Engbert
% Please cite: Collins, J.J., De Luca, C.J. Open-loop and closed-loop 
% control of posture: a random-walk analysis of center-of-pressure 
% trajectories. Exp. Brain. Res. 1993; 95:308-318
% Engbert, R. & Kliegl, R. Microsaccades keep the eyes' balance during 
% fixation. Psychological Science. 2004; 15:431-436.
%---------------------------------------------------------------------

N = length(x(:,1));
maxlag = round(N/4);
x1 = [x;x];
x2 = x;
r = zeros(1,maxlag);
for lag = 1:maxlag
  f=(lag:lag+N-1); 
  x11 = x1(lag:(lag+N-1),:);
  d = x11 - x2;
  r(lag) = mean(d(:,1).^2+d(:,2).^2);
end
lag = 1:maxlag;
rv = [lag; r];
