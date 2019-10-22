function a = tocmplx(a)
%TOCMPLX      Type cast to complex interval
%
%   z = tocmplx(a)
%

% written  10/16/98     S.M. Rump
% modified 11/30/98     S.M. Rump  improved speed
% modified 09/02/00     S.M. Rump  same as cintval
%                                  rounding unchanged after use
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
% modified 12/03/05     S.M. Rump  obsolete
%

  error('please replace "tocmplx" by "cintval"')
  