function [h0,h1] = hank103_jgh1(z)
%
%     

  h0 = complex(0);
  h1 = complex(0);
  ifexpon = 1;
  mex_id_ = 'hank103(i dcomplex[x], io dcomplex[x], io dcomplex[x], i int64_t[x])';
[h0, h1] = hank103_jgh(mex_id_, z, h0, h1, ifexpon, 1, 1, 1, 1);
end
%
%
%------------------------------------------------------
