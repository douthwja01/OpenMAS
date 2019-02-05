function y = test_h(x,index)
%TEST_H       Some test functions collected from Coconut
%
%   y = test_h(x,index)
%
%Input index specifies test function. Copyright see below.
%

% written  04/04/04     S.M. Rump
%

switch index
  
  case 1     % source see bottom of file
    y = x(3)-1 + x(1).^2 + x(2).^2 + (x(3)+x(4)).^2 + sin(x(3)).^2 + x(1).^2*x(2).^2 + ...
        x(4)-3 + sin(x(3)).^2 + (x(4)-1).^2 + x(2).^4 + x(3).^4 + (x(4)+x(1)).^2 + ...
        (x(1)-4 + sin(x(4)).^2 + x(2).^2*x(3).^2).^2 + sin(x(4)).^4;
      
  case 2     % source see bottom of file
    N = length(x);      % model problem: N = 1000, initial x=ones(N,1);
    I = 1:N-4;
    y = sum( (-4*x(I)+3.0).^2 ) + sum( ( x(I).^2 + 2*x(I+1).^2 + ...
              3*x(I+2).^2 + 4*x(I+3).^2 + 5*x(N).^2 ).^2 );
    
end


% function 1 taken from   http://www.sor.princeton.edu/~rvdb/ampl/nlmodels/cute/allinitu.mod
% # AMPL Model by Hande Y. Benson
% #
% # Copyright (C) 2001 Princeton University
% # All Rights Reserved
% #
% # Permission to use, copy, modify, and distribute this software and
% # its documentation for any purpose and without fee is hereby
% # granted, provided that the above copyright notice appear in all
% # copies and that the copyright notice and this
% # permission notice appear in all supporting documentation.                     
% 
% #   Source:
% #   N. Gould, private communication.
% 
% #   SIF input: Nick Gould, June 1990.
% 
% #   classification OUR2-AY-4-0
% 
% var x{1..4};
% 
% minimize f:
% x[3]-1 +
% x[1]^2+
% x[2]^2 + (x[3]+x[4])^2 +
% sin(x[3])^2 + x[1]^2*x[2]^2 + x[4]-3 +
% sin(x[3])^2 +
% (x[4]-1)^2 +
% (x[2]^2)^2+
% (x[3]^2 + (x[4]+x[1])^2)^2 +
% (x[1]-4 + sin(x[4])^2 + x[2]^2*x[3]^2)^2 +
% sin(x[4])^4;
% 
% solve;
% display f;
% display x;

% function 2 taken from   http://www.sor.princeton.edu/~rvdb/ampl/nlmodels/cute/bdqrtic.mod
% # AMPL Model by Hande Y. Benson
% #
% # Copyright (C) 2001 Princeton University
% # All Rights Reserved
% #
% # Permission to use, copy, modify, and distribute this software and
% # its documentation for any purpose and without fee is hereby
% # granted, provided that the above copyright notice appear in all
% # copies and that the copyright notice and this
% # permission notice appear in all supporting documentation.                     
% 
% #   Source: Problem 61 in
% #   A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
% #   "Performance of a multifrontal scheme for partially separable
% #   optimization",
% #   Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.
% 
% #   SIF input: Ph. Toint, Dec 1989.
% 
% #   classification SUR2-AN-V-0
% 
% param N:=1000;
% var x{1..N} := 1.0;
% 
% minimize f:
% sum {i in 1..N-4} (-4*x[i]+3.0)^2 + sum {i in 1..N-4} (x[i]^2+2*x[i+1]^2+3*x[i+2]^2+4*x[i+3]^2+5*x[N]^2)^2;
% 
% solve;
% display f;
% display x;