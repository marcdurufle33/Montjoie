function s = str2num_python(x)

  if (floor(x) == x)
    s = [num2str(x),'.0'];
  else
    s = num2str(x,12);
  end

