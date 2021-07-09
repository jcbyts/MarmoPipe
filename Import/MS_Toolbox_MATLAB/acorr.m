function [AC,L]=acorr(TS)
N=20;
fre =  2^(nextpow2(length(TS)) + 1);
F  =  fft(TS-mean(TS), fre);
AC  = ifft( F .* conj(F));
AC  =  real(AC(1:(N + 1))./AC(1));         % Retain non-negative lags.
L  =  [0:N]';