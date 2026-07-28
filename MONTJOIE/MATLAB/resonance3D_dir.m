function freq = resonance3D_dir(N,Nmax)

  freq = [];
  for n = 1:Nmax
    for m = 1:Nmax
      for k = 0:Nmax  
        freq = [freq, 0.5*sqrt(n*n+m*m+k*k)];
      end
    end
  end
  freq = sort(freq);
  freq = freq(1:N);
