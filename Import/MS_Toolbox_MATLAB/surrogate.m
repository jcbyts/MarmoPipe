function xsur = surrogate(x,SAMPLING) 
%---------------------------------------------------------------------
% Generation of surrogate data
% 2012, Ralf Engbert
% Please cite: Engbert, R., & Mergenthaler, K. (2006). Microsaccades 
% are triggered by low retinal image slip. Proceedings of the National 
% Academy of Sciences of the United States of America, 103, 7192-7197.
%---------------------------------------------------------------------

x0 = x(1,:);
v = vecvel(x,SAMPLING);
vsx = aaft(v(:,1))/SAMPLING;
vsy = aaft(v(:,2))/SAMPLING;
vsx(1) = vsx(1) + x0(1);
vsy(1) = vsy(1) + x0(2);
xs = cumsum(vsx);
ys = cumsum(vsy);
xsur = [xs' ys'];

