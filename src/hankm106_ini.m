function [ww,lw,ier] = hankm106_ini()
%
%     
  
  lw = 50000;
  ww = zeros(lw,1);
  ier= 0;

  mex_id_ = 'hank106datagen_r(io double[x], io int64_t[x], io int64_t[x])';
[ww, lw, ier] = hank106_jgh(mex_id_, ww, lw, ier, lw, 1, 1);

end
%
%
%------------------------------------------------------
%------------------------------------------------------
