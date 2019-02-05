void normalize(IntervalVector &x)
{
  IntervalVector n(x);

  for(int i = 0 ; i < x.size() ; i++)
  {
    Interval s = 0.;
    for(int j = 0 ; j < x.size() ; j++)
    {
      if(i == j) continue;
      s += sqr(x[j]);
    }

    Interval xi_pos = x[i] & Interval::POS_REALS;

    n[i] = 1. / sqrt((s/sqr(xi_pos)) + 1.);

    Interval xi_neg = x[i] & Interval::NEG_REALS;
    
    n[i] |= -1. / sqrt((s/sqr(xi_neg)) + 1.);
  }

  x = n;
}