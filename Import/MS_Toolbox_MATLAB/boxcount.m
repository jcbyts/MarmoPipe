function d = boxcount(xx,dx);
%-------------------------------------------------------------------
%  FUNCTION boxcount.m
%  Estimation of box-count;
%  Please cite: Engbert, R., & Mergenthaler, K. (2006) Microsaccades 
%  are triggered by low retinal image slip. Proceedings of the National 
%  Academy of Sciences of the United States of America, 103: 7192-7197.
%
%  (Version 1.0, April 2013)
%-------------------------------------------------------------------
%
%  INPUT: 
%
%  xx               epoch of trajectory
%  dx               edge length of boxes
%
%  OUTPUT:
%
%  d                number of boxes
%
%---------------------------------------------------------------------

y = xx(:,2);
x = xx(:,1);
x_min = min(x);
y_min = min(y);
x_max = max(x);
y_max = max(y);
MX = floor((x_max-x_min)/dx)+1;
MY = floor((y_max-y_min)/dx)+1;
box = zeros(MX,MY);
M = length(x);
for l=1:M
    i = floor( (x(l)-x_min)/dx ) + 1;
    j = floor( (y(l)-y_min)/dx ) + 1;
    box(i,j) = box(i,j) + 1;
end
d = length(find(box>0));

