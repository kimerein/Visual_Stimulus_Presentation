function index = binsearch(list, value)
% function index = binsearch(list, value)
% assumes list is sorted

idMin = 1;
idMax = length(list);

while idMin <= idMax
  idMid = floor(idMin + (idMax-idMin)/2);
  if value == list(idMid)
      index = idMid;
      return
  elseif value < list(idMid)
      idMax = idMid - 1;
  else
      idMin = idMid + 1;
  end
end
index = -1;