function err = getParameter(fichier, chaine)

  err = [];
  fid = fopen(fichier,'r');
  ligne = fgetl(fid);
  index = 1;
  while (feof(fid)==0)
    pos = strfind(ligne,chaine);
    if (length(pos) > 0)
      tab_numero = str2num(ligne(pos+length(chaine):end));
      err(index) = tab_numero;
      index = index+1;
    end
    ligne = fgetl(fid);
  end
  fclose(fid);
