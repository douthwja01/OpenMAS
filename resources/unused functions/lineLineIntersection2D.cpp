/*
   Determine the intersection point of two line segments
   Return FALSE if the lines don't intersect
*/
int LineIntersect(
double x1, double y1,
double x2, double y2,
double x3, double y3,
double x4, double y4,
double *x, double *y)
{
   double mua,mub;
   double denom,numera,numerb;

   denom  = (y4-y3) * (x2-x1) - (x4-x3) * (y2-y1);
   numera = (x4-x3) * (y1-y3) - (y4-y3) * (x1-x3);
   numerb = (x2-x1) * (y1-y3) - (y2-y1) * (x1-x3);

   /* Are the line coincident? */
   if (ABS(numera) < EPS && ABS(numerb) < EPS && ABS(denom) < EPS) {
      *x = (x1 + x2) / 2;
      *y = (y1 + y2) / 2;
      return(TRUE);
   }

   /* Are the line parallel */
   if (ABS(denom) < EPS) {
      *x = 0;
      *y = 0;
      return(FALSE);
   }

   /* Is the intersection along the the segments */
   mua = numera / denom;
   mub = numerb / denom;
   if (mua < 0 || mua > 1 || mub < 0 || mub > 1) {
      *x = 0;
      *y = 0;
      return(FALSE);
   }
   *x = x1 + mua * (x2 - x1);
   *y = y1 + mua * (y2 - y1);
   return(TRUE);
}