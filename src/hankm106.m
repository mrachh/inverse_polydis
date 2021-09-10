function [h0,h1] = hankm106(r,w,lw,intnum)
%
%     
  
  asize = size(r);
  tsize = numel(r);
  
  
  mex_id_ = 'hank106_wrap(i double[x], o dcomplex[x], o dcomplex[x], i double[x], i int64_t, i int64_t, i int64_t)';
[h0, h1] = hank106_jgh(mex_id_, r, w, lw, intnum, tsize, tsize, tsize, tsize, lw);

  h0 = reshape(h0,asize);
  h1 = reshape(h1,asize);

end
%
%
%------------------------------------------------------
