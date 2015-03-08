function [Sf, freq] = es_psd( sym_seq, dt )
% equal-symbol power spectral density
% reference: Voss, R.F., 1992, Evolution of Long-Rang Fractal Correlation
% and 1/f Noise in DNA Base Sequences
% dt is the time step (ms) for input data

N = length(sym_seq);
Fs = 1/dt*1000; % sampling frequency (Hz)


symbols = unique(sym_seq);
for i = 1:length(symbols)
    s = symbols(i);
    xdft = fft(sym_seq == s); % individual symbol-wise component
    
    % only need power estimates for the positive or negative frequencies.
    xdft = xdft(1:N/2+1);
    
    psdx = (1/(Fs*N)).*abs(xdft).^2;
    
    % In order to conserve the total power, multiply all frequencies that occur in both sets by a factor of 2.
    % Zero frequency (DC) and the Nyquist frequency do not occur twice.
    psdx(2:end-1) = 2*psdx(2:end-1);
    
    % accumulate for all the symbols
    if i == 1
        Sf = psdx;
    else
        Sf = Sf + psdx;
    end
end




freq = 0:Fs/length(sym_seq):Fs/2;


end