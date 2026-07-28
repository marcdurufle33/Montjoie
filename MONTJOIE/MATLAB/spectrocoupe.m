function []=spectrocoupe(y,nfft,fs,nwin,nover)

[B,f,t]=specgram(y,nfft,fs,nwin,nover);
%
% calculate amplitude 50dB down from maximum
%bmin=max(max(abs(B)))/300;
bmin=max(max(abs(B)))%/5000;
% plot top 50dB as image
%
imagesc(t,f,20*log10(abs(B)/bmin));



%bmin=2e-5
%imagesc(t,f,20*log10(abs(B)/bmin+eps))

% label plot
axis xy;
xlabel('Temps (s)');
ylabel('Frequence (Hz)');
%
% build and use a grey scale
%lgrays=zeros(100,3);
%for i=1:100
%    lgrays(i,:) = 1-i/100;
%end
%colormap(lgrays);
colormap(jet)
caxis([-90 0])
colorbar