function M=trace_fft(signal,F,col,posFigure)
% M=trace_fft(signal,F,col,posFigure)
% Displays the discrete Fourier transform of a signal, in dB, relatively to a maximal value.
% It will display it in the currently chosen figure.
% To have it in a new figure, 
% > figure, M=trace_fft(signal,F,col,posFigure)
% ----------------- INPUT VARIABLES -------------------
% signal : signal to process
% F : sampling frequency 
% col : (optional) color to use for the plot
% posFigure : (optional) position of the figure 
% ----------------- OUTPUT VARIABLES -------------------
% M : maximal value used for 0 dB

if nargin<3
    figure(gcf); hold on,
    col = 'b';
elseif nargin<4
     figure(gcf); hold on,
else
    h=figure(gcf);
   set(h,'Units','pixels')
    set(h,'Position',posFigure); 
end
% il y a un facteur sqrt(2 * pi)
toto=abs(fft(signal)/F/sqrt(2*pi));
M=max(toto);
%plot(linspace(-F/2,F/2,length(signal)),fftshift(20*log10(toto/M+eps)));
% normalisation 
plot(linspace(-F/2,F/2,length(signal)),fftshift(20*log10(toto/M+eps)),'Color',col);

xlim([0,F/2]);
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)')