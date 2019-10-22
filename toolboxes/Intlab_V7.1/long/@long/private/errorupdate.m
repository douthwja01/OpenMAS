function varargout = errorupdate(varargin)
%ERRORUPDATE  Addition of errors split by mantissa and exponent
%
%For internal use only

%
%Summation of errors for long numbers by mantissa and exponent and call
%
%  Err = errorupdate( f1,err1,e1 ,f2,err2,e2 ... )
%or
%  [ Errmant , Errexp ] = errorupdate( f1,err1,e1 ,f2,err2,e2 ... )
%
%where Err is structure with value Err.mant * beta^Err.exp and
%
%  Err >= sum( F_i * err_i * beta^e_i )
%
%and minimun one tripel, where
%  F_i = f_i    for nonnegative f_i (only sign of first component f_i(1) checked)
%  F_i = 1/f_i  otherwise.
%If err_i is structure, value is err_i.mant * beta^err_i.exp,
%otherwise err_i is the value itself.
%

% written  12/30/98     S.M. Rump
% modified 04/19/00     S.M. Rump  rounding unchanged after use
% modified 08/26/12     S.M. Rump  global variables removed
%

  INTLAB_LONG_BETA = getappdata(0,'INTLAB_LONG_BETA');

  setround(1)

  f = varargin{1};
  if f(1)<0
    f = 1./abs(f);
  end
  E = varargin{3};
  if isstruct(varargin{2})
    Err = varargin{2};
    Err.mant = f .* Err.mant;
    Err.exp = Err.exp + E;
  else
    Err.mant = f .* varargin{2};
    Err.exp = E;
  end

  arg = 3;
  while arg<nargin

    f = varargin{arg+1};
    if f(1)<0
      f = 1./abs(f);
    end
    e = varargin{arg+3};
    if isstruct(varargin{arg+2})
      err = varargin{arg+2};
      err.mant = f .* err.mant;
      err.exp = err.exp + e;
    else
      err.mant = f .* varargin{arg+2};
      err.exp = e;
    end
    arg = arg+3;

    % make sure sizes fit
    if ~isequal( length(Err.mant) , length(Err.exp) , ...
                 length(err.mant) , length(err.exp) )
      s = max([ length(Err.mant) length(Err.exp) ...
                length(err.mant) length(err.exp) ]);
      Err.mant = Err.mant .* ones(s,1);
      Err.exp = Err.exp .* ones(s,1);
      err.mant = err.mant .* ones(s,1);
      err.exp = err.exp .* ones(s,1);
    end

    % treat zero Err.mant
    indexE = ( Err.mant==0 );
    if any(indexE)
      Err.mant(indexE) = err.mant(indexE);
      Err.exp(indexE) = err.exp(indexE);
    end

    % indices of nonzero mantissas
    index = ~( indexE | ( err.mant==0 ) );

    % add errors for which Err.mant nonzero, Err.exp largest
    index1 = index & ( Err.exp>=err.exp );
    if any(index1)
      Err.mant(index1) = Err.mant(index1) + ...
        ( realmin + ...
          err.mant(index1) .* INTLAB_LONG_BETA.^(err.exp(index1)-Err.exp(index1)) );
    end

    % add errors for which Err.mant nonzero, err.exp largest
    index1 = index & ( Err.exp<err.exp );
    if any(index1)
      Err.mant(index1) = err.mant(index1) + ...
        ( realmin + ...
          Err.mant(index1) .* INTLAB_LONG_BETA.^(Err.exp(index1)-err.exp(index1)) );
      Err.exp(index1) = err.exp(index1);
    end

  end

  if nargout==1
    varargout{1} = Err;
  else
    varargout{1} = Err.mant;
    varargout{2} = Err.exp;
  end
