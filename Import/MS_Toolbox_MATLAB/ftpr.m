function xs = ftpr(x);
%-------------------------------------------------
%
% Compute Phase-Randomized surrogate data
% 
% 2012, Ralf Engbert
%
%-------------------------------------------------

N = length(x);
z = fftshift(fft(x));
phi = pi*(2*rand(1,N)-1);
z1 = z.*exp(i*phi);
z2re = real(z1+fliplr(z1))/2;
z2im = imag(z1-fliplr(z1))/2;
z2 = z2re + i*z2im;
xs = ifft(ifftshift(z2));
