function res = polynominit(param,see)
%POLANOMINIT  Initialization and defaults for polynom toolbox
%
%   polynominit(param,see)
%
%possible values for param:
%
%Default display of polynomials:
%  'DisplayUPolyVector'      Display univariate polynomials as vector, e.g.
%                               polynom p[x]  =
%                                    1     0    -2
%  'DisplayUPolySparse'      Sparse display of univariate polynomials, e.g.
%                               polynom p[x] =
%                                   2.0000 x^2
%                                  -1.0000
%  'DisplayUPoly'            res = 'DisplayUPolyVector'
%                                  'DisplayUPolySparse'
%
%Default evaluation of univariate interval polynomials (always matlab/polyval for no interval data):
%  'EvaluateUPolyHorner'     Evaluation by Horner's scheme, slow for large degree
%                               e.g. (x-3).*x+5
%  'EvaluateUPolyPower'      Evaluation by sums of scaled powers for degree>=20
%                               e.g. sum([1 -3 5].*[1 x x^2])
%  'EvaluateUPoly'           res = 'EvaluateUPolyHorner'
%                                  'EvaluateUPolyPower'
%
%Default evaluation of multivariate interval polynomials:
%  'EvaluateMPolyHorner'     Evaluation by recursive Horner's scheme
%  'EvaluateMPolyPower'      Evaluation by sums of scaled powers
%  'EvaluateMPoly'           res = 'EvaluateMPolyHorner'
%                                  'EvaluateMPolyPower'
%
%Default handling if polynomial does not depend on specified variable (for example, derivative,
%  scaling, shifting, transformation etc.)
%  'AccessVariableAuto'      Polynomial function executed without warning or error message even if polynomial 
%                               does not depend on specified variable
%  'AccessVariableWarn'      Polynomial function executed, but with warning
%  'AccessVariableError'     Polynomial function not executed but error message given (default)
%  'AccessVariable'          res = 'AccessVariableAuto'
%                                  'AccessVariableWarn'
%                                  'AccessVariableError'
%
%A corresponding message is printed unless it is explicitly suppressed with the (optional) parameter see=0.
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 08/26/12     S.M. Rump  global variables removed
% modified 11/13/12     S.M. Rump  reorganization
%

  if ( nargin==1 )
    see = 1;
  end
  
  switch lower(param)
    
    case 'init'
      setappdata(0,'INTLAB_POLYNOM_DISPLAY',0); % sparse display of univariate polynomials
      setappdata(0,'INTLAB_POLYNOM_UEVAL',1);   % Evaluation of univariate polynomials by Horner's scheme
      setappdata(0,'INTLAB_POLYNOM_MEVAL',1);   % Evaluation of multivariate polynomials by Horner's scheme
      setappdata(0,'INTLAB_POLYNOM_ACCESS_VARIABLE',1);  % Access to variable a polynomial does not depend on executed but with warning
      return
      
      
    case 'displayupoly'
      INTLAB_POLYNOM_DISPLAY = getappdata(0,'INTLAB_POLYNOM_DISPLAY');
      if INTLAB_POLYNOM_DISPLAY
        res = 'DisplayUPolyVector';
      else
        res = 'DisplayUPolySparse';
      end
      return
      
    case 'displayupolysparse'
      setappdata(0,'INTLAB_POLYNOM_DISPLAY',0);
      if see
        disp('===> Default sparse display of univariate polynomials')
      end
      return
      
    case 'displayupolyvector'
      setappdata(0,'INTLAB_POLYNOM_DISPLAY',1);
      if see
        disp('===> Default display of univariate polynomials as vector')
      end
      return
      
    case 'evaluateupoly'
      INTLAB_POLYNOM_UEVAL = getappdata(0,'INTLAB_POLYNOM_UEVAL');
      if INTLAB_POLYNOM_UEVAL
        res = 'EvaluateUPolyHorner';
      else
        res = 'EvaluateUPolyPower';
      end
      return
      
    case 'evaluateupolypower'
      setappdata(0,'INTLAB_POLYNOM_UEVAL',0);
      if see
        disp('===> Default evaluation univariate polynomials by sums of scaled powers')
      end
      return
      
    case 'evaluateupolyhorner'
      setappdata(0,'INTLAB_POLYNOM_UEVAL',1);
      if see
        disp('===> Default evaluation univariate polynomials by Horner''s scheme')
      end
      return
      
    case 'evaluatempoly'
      INTLAB_POLYNOM_MEVAL = getappdata(0,'INTLAB_POLYNOM_MEVAL');
      if INTLAB_POLYNOM_MEVAL
        res = 'EvaluateMPolyHorner';
      else
        res = 'EvaluateMPolyPower';
      end
      return
      
    case 'evaluatempolypower'
      setappdata(0,'INTLAB_POLYNOM_MEVAL',0);
      if see
        disp('===> Default evaluation multivariate polynomials by sums of scaled powers')
      end
      return
      
    case 'evaluatempolyhorner'
      setappdata(0,'INTLAB_POLYNOM_MEVAL',1);
      if see
        disp('===> Default evaluation multivariate polynomials by Horner''s scheme')
      end
      return
      
    case 'accessvariable'
      INTLAB_POLYNOM_ACCESS_VARIABLE = getappdata(0,'INTLAB_POLYNOM_ACCESS_VARIABLE');
      switch INTLAB_POLYNOM_ACCESS_VARIABLE
        case 0, res = 'AccessVariableAuto';
        case 1, res = 'AccessVariableWarn';
        case 2, res = 'AccessVariableError';
      end
      return
      
    case 'accessvariableauto'
      setappdata(0,'INTLAB_POLYNOM_ACCESS_VARIABLE',0);
      if see
        disp('===> Polynomial function executed without warning even if input polynomial does not depend on specified variable')
      end
      return
      
    case 'accessvariablewarn'
      setappdata(0,'INTLAB_POLYNOM_ACCESS_VARIABLE',1);
      if see
        disp('===> Polynomial function executed with warning if input polynomial does not depend on specified variable')
      end
      return
      
    case 'accessvariableerror'
      setappdata(0,'INTLAB_POLYNOM_ACCESS_VARIABLE',2);
      if see
        disp('===> Polynomial function: error message if input polynomial does not depend on specified variable')
      end
      return
      
      
    otherwise
      error('polynominit called with invalid argument')
      
  end
