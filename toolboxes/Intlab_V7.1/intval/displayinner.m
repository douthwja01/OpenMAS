function s = displayinner(yinf,ysup)
%DISPLAY      Command window display of inner inclusions (rigorous)
%
%If yinf,ysup is an inner inclusion computed by 
%     [x,yinf,ysup] = verifylss(A,b), 
%then for each 1<=i<=n the following is true:
%  There exist A1,A2 in A and b1,b2 in b such that
%(*)  (A1^-1 b1)_i <= yinf_i   and   ysup_i <= (A2^-1 b2)_i .
%
%To ensure rigor, an upper bound of yinf and a lower bound of ysup has to
%be displayed. This is done by
%     displayinner(yinf,ysup),
%where always infsup-notation is used. The output can be stored in a string by
%     str = displayinner(yinf,ysup)
%As for infsup display, the string shows always a row vector.
%
%The display by disp_ is not allowed for inner inclusions.
%
%Note that (*) is still true for yinf_i > ysup_i. In that case negative
%radii are displayed when using midrad format; thus infsup format is
%recommended for that case.
%

% written  06/15/11     S.M. Rump
% modified 06/15/11     S.M. Rump  string output, only infsup
%

  if isreal(yinf) & isreal(ysup)
    y = intval(yinf,ysup,'infsup');   % construct interval no matter of yinf<=ysup
    if nargout
      s = infsup(y,'inner');
    else
      infsup(y,'inner')
    end
  else
    error('Display of inner inclusions only for real arguments.')
  end
