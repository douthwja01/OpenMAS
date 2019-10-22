function plotintval(a,kmax)
%PLOTINTVAL   Plots real or complex interval
%
%  plotintval(a,kmax)
%
%Plots real 2xn interval vector as n boxes,
%  complex intervals are plotted as circles into Gaussian plane.
%Plots into current axis if hold on is active. Aspect ratio is set equal so that
%  squares and circles appear as such.
%Intervals of small diameter are plotted as asterisk.
%Parameter kmax for number of meshpoints for complex circles is optional, default is kmax=100.
%
%Example displaying wrapping effect:
%  a = midrad(0.8+0.7i,0.1); b = a; plotintval(a); hold on
%  for i=1:5, b = b*a, plotintval(b), end, hold off
%Another way to obtain the picture is
%  a = midrad(0.8+0.7i,0.1); for i=2:5, a(i) = a(i-1)*a(1); end, plotintval(a)
%A nice picture demonstrating overestimation of complex multiplication is generated as follows (data taken
%from Markus Neher: On the complex mean value form, SCAN 2002, Paris):
%
% kmax = 40; i = sqrt(-1); a=midrad(i,1); b=midrad(2+i,2); 
% plotintval(a*b); hold on
% phi = linspace(0,2*pi,kmax);
% [A,B] = meshgrid( mid(a)+rad(a)*exp(i*phi) , mid(b)+rad(b)*exp(i*phi) );
% plot(A.*B)
% hold off
%
%For another nice picture try a=midrad(1,1) and b=midrad(-1+i,1), presented by 
%Markus Grimmer from Wuppertal at the Dagstuhl meeting on Numerical Software with
%Result Verification, January 2003.
%

% written  08/06/02     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
% modified 11/20/05     S.M. Rump  fast check for rounding to nearest
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  if a.complex
    if nargin==1
      kmax = 100;
    end
    phi = linspace(0,2*pi,kmax)';
    mid = ones(kmax,1)*( a.mid(:).' );
    x = real(mid) + cos(phi)*a.rad(:)';
    y = imag(mid) + sin(phi)*a.rad(:)';
    index = ( a.rad(:)==0 );
    if any(index)
      plot( x,y, real(a.mid(index)),imag(a.mid(index)),'*' );
    else
      plot(x,y);
    end
    if ~ishold
      xmin = min(x(:));
      xmax = max(x(:));
      ymin = min(y(:));
      ymax = max(y(:));
      dx = 0.2*(xmax-xmin);
      dy = 0.2*(ymax-ymin);
      axis([xmin-dx xmax+dx ymin-dy ymax+dy]);
    end
  else
    ainf = a.inf;
    asup = a.sup;
    amin = min(ainf(:));
    amax = max(asup(:));
    d = 0.2*(amax-amin);
    if ( ndims(ainf)==2 ) & ( size(ainf,1)==2 )
      index = ( (asup(1,:)-ainf(1,:))<.1*d ) & ( (asup(2,:)-ainf(2,:))<.1*d );
      if any(index)      
        plot([ainf(1,:);asup(1,:);asup(1,:);ainf(1,:);ainf(1,:)],[ainf(2,:);ainf(2,:);asup(2,:);asup(2,:);ainf(2,:)], ainf(index),asup(index),'*' );
      else
        plot([ainf(1,:);asup(1,:);asup(1,:);ainf(1,:);ainf(1,:)],[ainf(2,:);ainf(2,:);asup(2,:);asup(2,:);ainf(2,:)]);
      end
      if ~ishold
        axis([amin-d amax+d amin-d amax+d]);
      end 
    else
      ainf = ainf(:);
      asup = asup(:);
      n = length(ainf);
      index = find( (asup-ainf)<.1*d );
      if any(index)      
        plot( [ainf asup]',[ 1:n ; 1:n ], [ ainf(index) asup(index)]',[index index]','*' );
      else
        plot([ainf asup]',[ 1:n ; 1:n ]);
      end
      if ~ishold
        axis([amin-d amax+d -1 n+1]);
      end 
    end   
  end
  axis equal
  
  setround(rndold)
