%Help file for INTLAB Version 3
%
%Additions/changes:
%
%  - Rigorous complex standard functions.
%  - Rigorous real standard functions revised, faster and more accurate.
%  - Suitable for Matlab student version with array length <=16384
%  - Fast argument reduction for |x|<=1e8.
%  - Built-in switch of rounding mode for PCs and Matlab Version 6 used.
%  - Verified inclusions for clustered eigenvalues and
%      invariant subspaces.
%  - Some auxiliary functions for long numbers.
%  - INTLAB input/output for multi-dimensional arrays.
%  - Standard functions for sparse input produce sparse output
%      iff f(0)=0 [otherwise result is full anyway, like sparse+1].
%  - Work directory defined in startup file.
%  - Function "isempty" for intervals redefined.
%  - For safety reasons, the data files storing local information
%      for and standard functions are appended with Matlab version
%      number. This is to make sure that always correct tables are
%      in use.
%  - Initialization procedure for standard functions much faster, about
%      30 minutes on a 300 Mhz PC.
%
%  Others:
%
%  - Global switch to control behaviour when real input interval
%      to standard function is out of range (see intvalinit)
%  - Type cast to sparse interval matrices added
%  - Type cast to full interval matrices added
%  - Functions length and size for \long added
%  - Some changes in subdirectory structure
%
%
%
%Further changes:
%
%  - Fix of one or the other bug, thanks to many users
%  - INTLAB logo
%
