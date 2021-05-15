namespace RobustPredicates
{
    internal static class ArithmeticFunctionsHelpers
    {
        internal static double Estimate(double[] e)
        {
            int eindex;

            double Q = e[0];
            for (eindex = 1; eindex < e.Length; eindex++)
            {
                Q += e[eindex];
            }
            return Q;
        }


        internal static int ScaleExpansionZeroeLim(double[] e, double b, double[] h)
        {
            int elen = e.Length;
            MacrosHelpers.Split(b, out double bhi, out double blo);
            MacrosHelpers.TwoProductPresplit(e[0], b, bhi, blo, out double Q, out double hh);
            int hindex = 0;
            if (hh != 0)
            {
                h[hindex++] = hh;
            }
            for (int eindex = 1; eindex < elen; eindex++)
            {
                double enow = e[eindex];
                MacrosHelpers.TwoProductPresplit(enow, b, bhi, blo, out double product1, out double product0);
                MacrosHelpers.TwoSum(Q, product0, out double sum, out hh);
                if (hh != 0)
                {
                    h[hindex++] = hh;
                }
                MacrosHelpers.FastTwoSum(product1, sum, out Q, out hh);
                if (hh != 0)
                {
                    h[hindex++] = hh;
                }
            }
            if ((Q != 0.0) || (hindex == 0))
            {
                h[hindex++] = Q;
            }
            return hindex;
        }

        internal static int FastExpansionSumZeroeLim(int elen, double[] e, int flen, double[] f, double[] h)
        {
            int findex;
            double Q;
            double Qnew;
            double hh;
            double enow = e[0];
            double fnow = f[0];
            int eindex = findex = 0;
            if ((fnow > enow) == (fnow > -enow))
            {
                Q = enow;
                enow = e[++eindex];
            }
            else
            {
                Q = fnow;
                fnow = f[++findex];
            }
            int hindex = 0;
            if ((eindex < elen) && (findex < flen))
            {
                if ((fnow > enow) == (fnow > -enow))
                {
                    MacrosHelpers.FastTwoSum(enow, Q, out Qnew, out hh);
                    enow = e[++eindex];
                }
                else
                {
                    MacrosHelpers.FastTwoSum(fnow, Q, out Qnew, out hh);
                    fnow = f[++findex];
                }
                Q = Qnew;
                if (hh != 0.0)
                {
                    h[hindex++] = hh;
                }
                while ((eindex < elen) && (findex < flen))
                {
                    if ((fnow > enow) == (fnow > -enow))
                    {
                        MacrosHelpers.TwoSum(Q, enow, out Qnew, out hh);
                        enow = e[++eindex];
                    }
                    else
                    {
                        MacrosHelpers.TwoSum(Q, fnow, out Qnew, out hh);
                        fnow = f[++findex];
                    }
                    Q = Qnew;
                    if (hh != 0.0)
                    {
                        h[hindex++] = hh;
                    }
                }
            }
            while (eindex < elen)
            {
                MacrosHelpers.TwoSum(Q, enow, out Qnew, out hh);
                enow = e[++eindex];
                Q = Qnew;
                if (hh != 0.0)
                {
                    h[hindex++] = hh;
                }
            }
            while (findex < flen)
            {
                MacrosHelpers.TwoSum(Q, fnow, out Qnew, out hh);
                fnow = f[++findex];
                Q = Qnew;
                if (hh != 0.0)
                {
                    h[hindex++] = hh;
                }
            }
            if ((Q != 0.0) || (hindex == 0))
            {
                h[hindex++] = Q;
            }

            return hindex;
        }
    }
}
