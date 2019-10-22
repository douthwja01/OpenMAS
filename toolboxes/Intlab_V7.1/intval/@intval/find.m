function [I,J,V] = find(a)
%FIND         Implements  find(a)  for sparse interval matrix
%
%   I = find(a)
%   [I,J] = find(a)
%   [I,J,V] = find(a)
%
%Functionality as in Matlab.
%

% written  08/09/02     S.M. Rump 
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    Matlab bug for empty input fixed
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 09/28/08     S.M. Rump  check for rounding to nearest improved
% modified 10/03/12     S.M. Rump  spones for ND-arrays
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if nargout<=1
    if a.complex
      if isequal(a.rad,0)
        I = find(a.mid);
      else
        I = find(abs(a.mid)+abs(a.rad));
      end
    else
      I = find(abs(a.inf)+abs(a.sup));
    end
  elseif nargout==2
    if a.complex
      if isequal(a.rad,0)
        [I,J] = find(a.mid);
      else
        [I,J] = find(abs(a.mid)+abs(a.rad));
      end
    else
      [I,J] = find(abs(a.inf)+abs(a.sup));
    end
  elseif nargout==3  
    if a.complex
      if prod(size(a.mid))<2^31
        if isequal(a.rad,0)
          [I,J] = find(a.mid);
        else
          [I,J] = find(abs(a.mid)+abs(a.rad));
        end
        if isempty(I)                   % cures Matlab bug
          K = I;
        elseif isempty(J)               % cures Matlab bug
          K = J;
        else
          K = sub2ind(size(a.mid),I,J);
        end
        if isequal(a.rad,0)
          V = intval(full(a.mid(K)),0,'midrad');
        else
          V = intval(full(a.mid(K)),full(a.rad(K)),'midrad');
        end
      else
        smid = spones(a.mid(:,:));
        srad = spones(a.rad(:,:));
        if ~isequal(a.rad,0) & ~isequal(smid,srad)  % take care of <0,x> and <x,0> components
          smidrad = spones(a.mid(:,:)).*spones(a.rad(:,:));
          [I J] = find(spones(a.mid(:,:))+spones(a.rad(:,:)));
          zmid = srad-smidrad;
          zrad = smid-smidrad;
          a.mid = a.mid + zmid;               % set <0,x> components to <1,x>
          a.rad = a.rad + zrad;               % set <x,0> components to <x,1>
          a.mid = nonzeros(a.mid) + ( 1 - nonzeros(spones(a.mid(:,:))+zmid) );
          a.rad = nonzeros(a.rad) + ( 1 - nonzeros(spones(a.rad(:,:))+zrad) );
          V = intval(a.mid,a.rad,'midrad');
        else
          [I J] = find(smid);
          if isequal(a.rad,0)
            V = intval(full(nonzeros(a.mid)),0,'midrad');
          else
            V = intval(full(nonzeros(a.mid)),full(nonzeros(a.rad)),'midrad');
          end
        end
      end
    else
      if prod(size(a.inf))<2^31
        [I,J] = find(abs(a.inf)+abs(a.sup));
        if isempty(I)                     % cures Matlab bug
          K = I;
        elseif isempty(J)                 % cures Matlab bug
          K = J;
        else
          K = sub2ind(size(a.inf),I,J);
        end
        V = intval(a.inf(K),a.sup(K),'infsup');
      else
        sinf = spones(a.inf(:,:));
        ssup = spones(a.sup(:,:));
        if ~isequal(sinf,ssup)                % take care of [0,x] and [x,0] components
          sinfsup = spones(a.inf(:,:)).*spones(a.sup(:,:));
          [I J] = find(spones(a.inf(:,:))+spones(a.sup(:,:)));
          zinf = ssup-sinfsup;
          zsup = sinf-sinfsup;
          a.inf = a.inf + zinf;               % set [0,x] components to [1,x]
          a.sup = a.sup + zsup;               % set [x,0] components to [x,1]
          a.inf = nonzeros(a.inf) + ( 1 - nonzeros(spones(a.inf(:,:))+zinf) );
          a.sup = nonzeros(a.sup) + ( 1 - nonzeros(spones(a.sup(:,:))+zsup) );
          V = intval(a.inf,a.sup,'infsup');
        else
          [I J] = find(sinf);
          V = intval(nonzeros(a.inf),nonzeros(a.sup),'infsup');
        end
      end
    end
  end
  
  if rndold
    setround(rndold)
  end
