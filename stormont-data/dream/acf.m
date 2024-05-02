function [acf] = acf(y,numLags)

nFFT = 2^(nextpow2(length(y))+1);
F = fft(y-mean(y),nFFT);
F = F.*conj(F);
acf = ifft(F);
acf = acf(1:(numLags+1)); % Retain non-negative lags
acf = acf./acf(1); % Normalize
%%acf = acf./length(y);
acf = real(acf(2:end)); % don't return peak of 1


