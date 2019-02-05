%Help file for INTLAB Version 2
%
%Changes:
%
%  Input/output:
%
%  - INTLAB output is rigorous [ Note that INTLAB V1 input is already rigorous
%      when using character constants (e.g. intval('.1')) ]
%  - rigorous standard functions
%  - Interval display by uncertainty, e.g. 3.14159_
%  - Interval input by string, also with tolerances (e.g. '3.14_',
%      '[3,4]' or '<-3.1e2,0.01_>' or '<3-4i,.1> etc.
%  - Permanent switch of display mode of intervals thru intvalinit
%  - Function initvar for gradients replaced by gradientinit (replace
%      by global change, exactly the same functionality); was necessary
%      because of ambiguouity with initialization of slope variables
%
%
%  Standard functions:
%
%  - Rigorous standard functions incorporated
%
%
%  Slopes:
%
%  - Slope toolbox added
%
%
%  Long precision:
%
%  - Long toolbox added
%
%
%  Others:
%
%  - Simplified call of verifynlss for one-dimensional nonlinear functions
%  - Extended arithmetic including +/-infinity, e.g. division by
%      zero intervals etc.
%  - Some changes in subdirectory structure
%
%
%
%Further changes:
%
%  - Linear system solver in 2 stages with sharp interval hull for
%      preconditioned system
%  - Components of arrays of intervals, gradients... may be erased by
%      assignment of an empty set, e.g. A=intval(ones(5)); A(2,:)=[];
%  - Functions SetRoundDown, SetRoundUp and SetRoundNear replaced by
%      one routine setround(rnd)
%  - Functions "in", "in0", "isnan" deliver array of 0/1 (like Matlab)
%  - Functions pred/succ with second parameter for k-th predecessor/successor
%  - Function power10tobinary with increased exponent range, changes also
%      power10tobinary.mat
%  - Fix of one or the other bug, thanks to many users
%
%New functions despite the above:
%
%  - mig          Mignitude for points and intervals
%  - compmat      Ostrowski's comparison matrix
%  - band         Extract band out of matrix
%  - bandwidth    Get lower and upper bandwidth
%  - midradcmplx  Always produces complex interval
%  - binom        Binomial coefficient
%  - distbits     Distance in bits of two reals
%  - getround     Get current rounding mode
%  - gregk...     Some sparse test matrices from Gregory/Karney test set
%  - perms_       Fast version of Matlab perms
%  - randint      Random integers in specified range
%  - random       Random numbers uniformly distrubuted in [-1,1] or
%                   within specified range
%  - randomc      Complex random numbers with real and imaginary part as random
%  - relerr       Relative error for scalars, vectors and matrices
%
