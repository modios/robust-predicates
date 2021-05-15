using System;

namespace RobustPredicates
{
    public static class Orient2D
    {
        private static double Orient2dAdapt(double[] pa, double[] pb, double[] pc, double detsum)
        {
            double acx = pa[0] - pc[0];
            double bcx = pb[0] - pc[0];
            double acy = pa[1] - pc[1];
            double bcy = pb[1] - pc[1];

            MacrosHelpers.TwoProduct(acx, bcy, out double detleft, out double detlefttail);
            MacrosHelpers.TwoProduct(acy, bcx, out double detright, out double detrighttail);
            double[] B = new double[4];
            MacrosHelpers.TwoTwoDiff(detleft, detlefttail, detright, detrighttail, out B[3], out B[2], out B[1], out B[0]);

            double det = ArithmeticFunctionsHelpers.Estimate(B);
            double errbound = MacrosHelpers.CcwerrboundB * detsum;
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            MacrosHelpers.TwoDiffTail(pa[0], pc[0], acx, out double acxtail);
            MacrosHelpers.TwoDiffTail(pb[0], pc[0], bcx, out double bcxtail);
            MacrosHelpers.TwoDiffTail(pa[1], pc[1], acy, out double acytail);
            MacrosHelpers.TwoDiffTail(pb[1], pc[1], bcy, out double bcytail);

            if ((acxtail == 0.0) && (acytail == 0.0)
                && (bcxtail == 0.0) && (bcytail == 0.0))
            {
                return det;
            }

            errbound = MacrosHelpers.CcwerrboundC * detsum + MacrosHelpers.Resulterrbound * Math.Abs(det);
            det += (acx * bcytail + bcy * acxtail)
                 - (acy * bcxtail + bcx * acytail);
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            double[] u = new double[4];
            MacrosHelpers.TwoProduct(acxtail, bcy, out double s1, out double s0);
            MacrosHelpers.TwoProduct(acytail, bcx, out double t1, out double t0);
            MacrosHelpers.TwoTwoDiff(s1, s0, t1, t0, out u[3], out u[2], out u[1], out u[0]);

            double[] C1 = new double[8];
            int C1length = ArithmeticFunctionsHelpers.FastExpansionSumZeroeLim(4, B, 4, u, C1);

            MacrosHelpers.TwoProduct(acx, bcytail, out s1, out s0);
            MacrosHelpers.TwoProduct(acy, bcxtail, out t1, out t0);
            MacrosHelpers.TwoTwoDiff(s1, s0, t1, t0, out u[3], out u[2], out u[1], out u[0]);
            double[] C2 = new double[12];
            int C2length = ArithmeticFunctionsHelpers.FastExpansionSumZeroeLim(C1length, C1, 4, u, C2);

            MacrosHelpers.TwoProduct(acxtail, bcytail, out s1, out s0);
            MacrosHelpers.TwoProduct(acytail, bcxtail, out t1, out t0);
            MacrosHelpers.TwoTwoDiff(s1, s0, t1, t0, out u[3], out u[2], out u[1], out u[0]);

            double[] D = new double[16];
            int Dlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroeLim(C2length, C2, 4, u, D);

            return D[Dlength - 1];
        }

        public static double Fast(double[] pa, double[] pb, double[] pc)
        {
            double acx = pa[0] - pc[0];
            double bcx = pb[0] - pc[0];
            double acy = pa[1] - pc[1];
            double bcy = pb[1] - pc[1];
            return acx * bcy - acy * bcx;
        }

        internal static double Exact(double[] pa, double[] pb, double[] pc)
        {
            MacrosHelpers.TwoProduct(pa[0], pb[1], out double axby1, out double axby0);
            MacrosHelpers.TwoProduct(pa[0], pc[1], out double axcy1, out double axcy0);
            double[] aterms = new double[4];
            MacrosHelpers.TwoTwoDiff(axby1, axby0, axcy1, axcy0, out aterms[3], out aterms[2], out aterms[1], out aterms[0]);

            MacrosHelpers.TwoProduct(pb[0], pc[1], out double bxcy1, out double bxcy0);
            MacrosHelpers.TwoProduct(pb[0], pa[1], out double bxay1, out double bxay0);

            double[] bterms = new double[4];
            MacrosHelpers.TwoTwoDiff(bxcy1, bxcy0, bxay1, bxay0, out bterms[3], out bterms[2], out bterms[1], out bterms[0]);

            MacrosHelpers.TwoProduct(pc[0], pa[1], out double cxay1, out double cxay0);
            MacrosHelpers.TwoProduct(pc[0], pb[1], out double cxby1, out double cxby0);


            double[] cterms = new double[4];
            MacrosHelpers.TwoTwoDiff(cxay1, cxay0, cxby1, cxby0, out cterms[3], out cterms[2], out cterms[1], out cterms[0]);

            double[] v = new double[8];
            double[] w = new double[12];
            int vlen = ArithmeticFunctionsHelpers.FastExpansionSumZeroeLim(4, aterms, 4, bterms, v);
            int wlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroeLim(vlen, v, 4, cterms, w);

            return w[wlength - 1];
        }

        internal static double Slow(double[] pa, double[] pb, double[] pc)
        {
            MacrosHelpers.TwoDiff(pa[0], pc[0], out double acx, out double acxtail);
            MacrosHelpers.TwoDiff(pa[1], pc[1], out double acy, out double acytail);
            MacrosHelpers.TwoDiff(pb[0], pc[0], out double bcx, out double bcxtail);
            MacrosHelpers.TwoDiff(pb[1], pc[1], out double bcy, out double bcytail);
            double[] axby = new double[8];

            MacrosHelpers.TwoTwoProduct(acx, acxtail, bcy, bcytail,
                            out axby[7], out axby[6], out axby[5], out axby[4],
                            out axby[3], out axby[2], out axby[1], out axby[0]);
            double negate = -acy;
            double negatetail = -acytail;
            double[] bxay = new double[8];
            MacrosHelpers.TwoTwoProduct(bcx, bcxtail, negate, negatetail,
                            out bxay[7], out bxay[6], out bxay[5], out bxay[4],
                            out bxay[3], out bxay[2], out bxay[1], out bxay[0]);

            double[] deter = new double[16];
            int deterlen = ArithmeticFunctionsHelpers.FastExpansionSumZeroeLim(8, axby, 8, bxay, deter);

            return deter[deterlen - 1];
        }

        public static double Robust(double[] pa, double[] pb, double[] pc)
        {

            double detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
            double detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
            double det = detleft - detright;
            double detsum;

            if (detleft > 0.0)
            {
                if (detright <= 0.0)
                {
                    return det;
                }
                else
                {
                    detsum = detleft + detright;
                }
            }
            else if (detleft < 0.0)
            {
                if (detright >= 0.0)
                {
                    return det;
                }
                else
                {
                    detsum = -detleft - detright;
                }
            }
            else
            {
                return det;
            }

            double errbound = MacrosHelpers.CcwerrboundA * detsum;
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            return Orient2dAdapt(pa, pb, pc, detsum);
        }
    }
}
