function freq = resonance(N,Nmax)

  freq = [];
  for n = 0:Nmax
    for m = 0:Nmax
      freq = [freq, 0.5*sqrt(n*n+m*m)];
    end
  end
  freq = sort(freq);
  freq = freq(1:N);
