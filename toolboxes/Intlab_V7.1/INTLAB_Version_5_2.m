%Help file for INTLAB Version 5.2
%
% - new algorithm AccSum for accurate summation: faithful, to nearest, K-fold and inclusion
% - new algorithm AccDot for accurate dot product: faithful, to nearest, K-fold and inclusion
% - option intvalinit('RealStdFctsExcptnIgnore') to ignore input arguments out of range (thanks to Arnold)
% - meaning of isempty/isempty_ for intervals changed: This was necessary
%     to ensure that the interval function 'isempty' is the natural
%     extension of 'isempty' for the point function. That means,
%     isempty(X) for an interval quantity checks whether X is [].
% - parameters of functions sum_ and dot_ adapted to appeared paper
% - improved performance for setround (thanks to Jˆrg Kubitz, Hannover)
% - improved performance for gradient/mtimes (thanks to Jˆrg Kubitz, Hannover)
% - improved performance for hessian/mtimes
% - intval(0)^0 changed to 1 (thanks to Jˆrg Kubitz, Hannover)
% - corrections in long toolbox (thanks to Nobito Yamamoto and Nozomu Matsuda, Tokyo)
% - intval/power corrected (thanks to John Pryce, Swindon)
% - note added to intval/rad (thanks to GÅEter Mayer, Rostock)
% - workpath omitted and directory unchanged in startintlab
% - functions ldivide added to all data types (thanks to Jˆrg Kubitz, Hannover)
% - files Version5.2 etc. renamed into Version5_2 to allow help Version5_2
% - tocmplx obsolete (replaced by cintval)
% - function find gradient, hessian and slope
% - functions any, all, logical for intval, gradient, hessian
%
