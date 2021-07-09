function y = smoothdata(x);
%--------------------------------------------------------------
% Smoothes data according to a 5-point moving average procedure
% 2012, Ralf Engbert
%--------------------------------------------------------------

y0 = x(1,:);
v = vecvel(x,1);
y = cumsum([y0; v]);
