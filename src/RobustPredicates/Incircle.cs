namespace RobustPredicates
{
    public static class InCirlce
    {
        public static double Fast(double[] pa, double[] pb, double[] pc, double[] pd)
        {
            double adx = pa[0] - pd[0];
            double ady = pa[1] - pd[1];
            double bdx = pb[0] - pd[0];
            double bdy = pb[1] - pd[1];
            double cdx = pc[0] - pd[0];
            double cdy = pc[1] - pd[1];

            double abdet = adx * bdy - bdx * ady;
            double bcdet = bdx * cdy - cdx * bdy;
            double cadet = cdx * ady - adx * cdy;
            double alift = adx * adx + ady * ady;
            double blift = bdx * bdx + bdy * bdy;
            double clift = cdx * cdx + cdy * cdy;

            return alift * bcdet + blift * cadet + clift * abdet;
        }


    }
}
