using System;

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

        internal static double Exact(double[] pa, double[] pb, double[] pc, double[] pd)
        {
            MacrosHelpers.TwoProduct(pa[0], pb[1], out double axby1, out double axby0);
            MacrosHelpers.TwoProduct(pb[0], pa[1], out double bxay1, out double bxay0);

            var ab = new double[4];
            MacrosHelpers.TwoTwoDiff(axby1, axby0, bxay1, bxay0, out ab[3], out ab[2], out ab[1], out ab[0]);

            var bc = new double[4];
            MacrosHelpers.TwoProduct(pb[0], pc[1], out double bxcy1, out double bxcy0);
            MacrosHelpers.TwoProduct(pc[0], pb[1], out double cxby1, out double cxby0);
            MacrosHelpers.TwoTwoDiff(bxcy1, bxcy0, cxby1, cxby0, out bc[3], out bc[2], out bc[1], out bc[0]);

            var cd = new double[4];
            MacrosHelpers.TwoProduct(pc[0], pd[1], out double cxdy1, out double cxdy0);
            MacrosHelpers.TwoProduct(pd[0], pc[1], out double dxcy1, out double dxcy0);
            MacrosHelpers.TwoTwoDiff(cxdy1, cxdy0, dxcy1, dxcy0, out cd[3], out cd[2], out cd[1], out cd[0]);

            var da = new double[4];
            MacrosHelpers.TwoProduct(pd[0], pa[1], out double dxay1, out double dxay0);
            MacrosHelpers.TwoProduct(pa[0], pd[1], out double axdy1, out double axdy0);
            MacrosHelpers.TwoTwoDiff(dxay1, dxay0, axdy1, axdy0, out da[3], out da[2], out da[1], out da[0]);

            var ac = new double[4];
            MacrosHelpers.TwoProduct(pa[0], pc[1], out double axcy1, out double axcy0);
            MacrosHelpers.TwoProduct(pc[0], pa[1], out double cxay1, out double cxay0);
            MacrosHelpers.TwoTwoDiff(axcy1, axcy0, cxay1, cxay0, out ac[3], out ac[2], out ac[1], out ac[0]);

            var bd = new double[4];
            MacrosHelpers.TwoProduct(pb[0], pd[1], out double bxdy1, out double bxdy0);
            MacrosHelpers.TwoProduct(pd[0], pb[1], out double dxby1, out double dxby0);
            MacrosHelpers.TwoTwoDiff(bxdy1, bxdy0, dxby1, dxby0, out bd[3], out bd[2], out bd[1], out bd[0]);

            var temp8 = new double[8];
            var templen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(4, cd, 4, da, temp8);
            var cda = new double[12];
            var cdalen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(templen, temp8, 4, ac, cda);
            var dab = new double[12];
            templen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(4, da, 4, ab, temp8);
            var dablen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(templen, temp8, 4, bd, dab);

            for (int i = 0; i < 4; i++)
            {
                bd[i] = -bd[i];
                ac[i] = -ac[i];
            }

            var abc = new double[12];
            templen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(4, ab, 4, bc, temp8);
            var abclen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(templen, temp8, 4, ac, abc);
            templen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(4, bc, 4, cd, temp8);

            var bcd = new double[12];
            var bcdlen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(templen, temp8, 4, bd, bcd);

            var det24x = new double[24];
            var xlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bcdlen, bcd, pa[0], det24x);
            double[] det48x = new double[48];
            xlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(xlen, det24x, pa[0], det48x);
            double[] det24y = new double[24];
            var ylen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bcdlen, bcd, pa[1], det24y);
            double[] det48y = new double[48];
            ylen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(ylen, det24y, pa[1], det48y);
            double[] adet = new double[96];
            var alen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(xlen, det48x, ylen, det48y, adet);

            xlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cdalen, cda, pb[0], det24x);
            xlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(xlen, det24x, -pb[0], det48x);
            ylen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cdalen, cda, pb[1], det24y);
            ylen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(ylen, det24y, -pb[1], det48y);
            double[] bdet = new double[96];
            var blen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(xlen, det48x, ylen, det48y, bdet);

            xlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(dablen, dab, pc[0], det24x);
            xlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(xlen, det24x, pc[0], det48x);
            ylen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(dablen, dab, pc[1], det24y);
            ylen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(ylen, det24y, pc[1], det48y);
            var cdet = new double[96];
            var clen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(xlen, det48x, ylen, det48y, cdet);

            xlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(abclen, abc, pd[0], det24x);
            xlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(xlen, det24x, -pd[0], det48x);
            ylen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(abclen, abc, pd[1], det24y);
            ylen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(ylen, det24y, -pd[1], det48y);
            double[] ddet = new double[96];
            var dlen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(xlen, det48x, ylen, det48y, ddet);

            double[] abdet = new double[192];
            var ablen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(alen, adet, blen, bdet, abdet);
            double[] cddet = new double[192];
            var cdlen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(clen, cdet, dlen, ddet, cddet);
            double[] deter = new double[384];
            var deterlen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(ablen, abdet, cdlen, cddet, deter);
            return deter[deterlen - 1];
        }

        private static double Adapt(double[] pa, double[] pb, double[] pc, double[] pd, double permanent)
        {
            var adx = pa[0] - pd[0];
            var bdx = pb[0] - pd[0];
            var cdx = pc[0] - pd[0];
            var ady = pa[1] - pd[1];
            var bdy = pb[1] - pd[1];
            var cdy = pc[1] - pd[1];

            MacrosHelpers.TwoProduct(bdx, cdy, out double bdxcdy1, out double bdxcdy0);
            MacrosHelpers.TwoProduct(cdx, bdy, out double cdxbdy1, out double cdxbdy0);
            var bc = new double[4];
            MacrosHelpers.TwoTwoDiff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, out bc[3], out bc[2], out bc[1], out bc[0]);

            double[] axbc = new double[8];
            var axbclen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, bc, adx, axbc);
            double[] axxbc = new double[16];
            var axxbclen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(axbclen, axbc, adx, axxbc);
            double[] aybc = new double[8];
            var aybclen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, bc, ady, aybc);
            double[] ayybc = new double[16];
            var ayybclen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(aybclen, aybc, ady, ayybc);
            double[] adet = new double[32];
            var alen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(axxbclen, axxbc, ayybclen, ayybc, adet);

            MacrosHelpers.TwoProduct(cdx, ady, out double cdxady1, out double cdxady0);
            MacrosHelpers.TwoProduct(adx, cdy, out double adxcdy1, out double adxcdy0);
            var ca = new double[4];
            MacrosHelpers.TwoTwoDiff(cdxady1, cdxady0, adxcdy1, adxcdy0, out ca[3], out ca[2], out ca[1], out ca[0]);

            double[] bxca = new double[8];
            var bxcalen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, ca, bdx, bxca);
            double[] bxxca = new double[16];
            var bxxcalen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bxcalen, bxca, bdx, bxxca);
            double[] byca = new double[8];
            var bycalen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, ca, bdy, byca);
            double[] byyca = new double[16];
            var byycalen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bycalen, byca, bdy, byyca);
            double[] bdet = new double[32];
            var blen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(bxxcalen, bxxca, byycalen, byyca, bdet);



            MacrosHelpers.TwoProduct(adx, bdy, out double adxbdy1, out double adxbdy0);
            MacrosHelpers.TwoProduct(bdx, ady, out double bdxady1, out double bdxady0);
            var ab = new double[4];
            MacrosHelpers.TwoTwoDiff(adxbdy1, adxbdy0, bdxady1, bdxady0, out ab[3], out ab[2], out ab[1], out ab[0]);

            double[] cxab = new double[8];
            var cxablen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, ab, cdx, cxab);
            double[] cxxab = new double[16];
            var cxxablen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cxablen, cxab, cdx, cxxab);

            double[] cyab = new double[8];
            var cyablen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, ab, cdy, cyab);
            double[] cyyab = new double[16];
            var cyyablen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cyablen, cyab, cdy, cyyab);
            double[] cdet = new double[32];
            var clen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(cxxablen, cxxab, cyyablen, cyyab, cdet);

            double[] abdet = new double[64];
            var ablen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(alen, adet, blen, bdet, abdet);
            double[] fin1 = new double[1152];
            var finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(ablen, abdet, clen, cdet, fin1);

            var det = ArithmeticFunctionsHelpers.Estimate(finlength, fin1);
            var errbound = MacrosHelpers.IccerrboundB * permanent;

            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            MacrosHelpers.TwoDiffTail(pa[0], pd[0], adx, out double adxtail);
            MacrosHelpers.TwoDiffTail(pa[1], pd[1], ady, out double adytail);
            MacrosHelpers.TwoDiffTail(pb[0], pd[0], bdx, out double bdxtail);
            MacrosHelpers.TwoDiffTail(pb[1], pd[1], bdy, out double bdytail);
            MacrosHelpers.TwoDiffTail(pc[0], pd[0], cdx, out double cdxtail);
            MacrosHelpers.TwoDiffTail(pc[1], pd[1], cdy, out double cdytail);

            if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
                && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0))
            {
                return det;
            }

            errbound = MacrosHelpers.IccerrboundC * permanent + MacrosHelpers.Resulterrbound * Math.Abs(det);
            det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail)
                                               - (bdy * cdxtail + cdx * bdytail))
                    + 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
                 + ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail)
                                               - (cdy * adxtail + adx * cdytail))
                    + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
                 + ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail)
                                               - (ady * bdxtail + bdx * adytail))
                    + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }


            var finnow = fin1;
            var fin2 = new double[1152];
            var finother = fin2;
            double[] aa = new double[4];
            double[] bb = new double[4];
            double[] cc = new double[4];

            if ((bdxtail != 0.0) || (bdytail != 0.0)
                || (cdxtail != 0.0) || (cdytail != 0.0))
            {
                MacrosHelpers.Square(adx, out double adxadx1, out double adxadx0);
                MacrosHelpers.Square(ady, out double adyady1, out double adyady0);
                MacrosHelpers.TwoTwoSum(adxadx1, adxadx0, adyady1, adyady0, out aa[3], out aa[2], out aa[1], out aa[0]);
            }
            if ((cdxtail != 0.0) || (cdytail != 0.0)
                || (adxtail != 0.0) || (adytail != 0.0))
            {
                MacrosHelpers.Square(bdx, out double bdxbdx1, out double bdxbdx0);
                MacrosHelpers.Square(bdy, out double bdybdy1, out double bdybdy0);
                MacrosHelpers.TwoTwoSum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, out bb[3], out bb[2], out bb[1], out bb[0]);
            }
            if ((adxtail != 0.0) || (adytail != 0.0)
                || (bdxtail != 0.0) || (bdytail != 0.0))
            {
                MacrosHelpers.Square(cdx, out double cdxcdx1, out double cdxcdx0);
                MacrosHelpers.Square(cdy, out double cdycdy1, out double cdycdy0);
                MacrosHelpers.TwoTwoSum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, out cc[3], out cc[2], out cc[1], out cc[0]);
            }


            double[] finswap;
            double[] temp16a = new double[16];
            double[] temp32a = new double[32];
            double[] temp16b = new double[16];
            double[] temp32b = new double[32];
            double[] temp16c = new double[16];
            double[] temp48 = new double[48];
            double[] axtbc = new double[8];
            int axtbclen = 0;
            if (adxtail != 0.0)
            {
                axtbclen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, bc, adxtail, axtbc);
                var temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(axtbclen, axtbc, 2.0 * adx, temp16a);
                double[] axtcc = new double[8];
                var axtcclen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, cc, adxtail, axtcc);
                var temp16blen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(axtcclen, axtcc, bdy, temp16b);
                double[] axtbb = new double[8];
                var axtbblen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, bb, adxtail, axtbb);
                var temp16clen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(axtbblen, axtbb, -cdy, temp16c);

                var temp32alen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                        temp16blen, temp16b, temp32a);

                var temp48len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16clen, temp16c,
                                                        temp32alen, temp32a, temp48);
                finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                        temp48, finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap;
            }

            int aytbclen = 0;
            double[] aytbc = new double[8];
            if (adytail != 0.0)
            {
                aytbclen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, bc, adytail, aytbc);
                var temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(aytbclen, aytbc, 2.0 * ady, temp16a);
                double[] aytbb = new double[8];
                var aytbblen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, bb, adytail, aytbb);
                var temp16blen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(aytbblen, aytbb, cdx, temp16b);
                double[] aytcc = new double[8];
                var aytcclen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, cc, adytail, aytcc);
                var temp16clen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(aytcclen, aytcc, -bdx, temp16c);
                var temp32alen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                        temp16blen, temp16b, temp32a);
                var temp48len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16clen, temp16c,
                                                        temp32alen, temp32a, temp48);
                finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                        temp48, finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap;
            }

            int bxtcalen = 0;
            double[] bxtca = new double[8];
            if (bdxtail != 0.0)
            {
                bxtcalen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, ca, bdxtail, bxtca);
                var temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bxtcalen, bxtca, 2.0 * bdx, temp16a);
                double[] bxtaa = new double[8];
                var bxtaalen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, aa, bdxtail, bxtaa);
                var temp16blen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bxtaalen, bxtaa, cdy, temp16b);
                double[] bxtcc = new double[8];
                var bxtcclen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, cc, bdxtail, bxtcc);
                var temp16clen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bxtcclen, bxtcc, -ady, temp16c);
                var temp32alen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                        temp16blen, temp16b, temp32a);
                var temp48len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16clen, temp16c,
                                                        temp32alen, temp32a, temp48);
                finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                        temp48, finother);
                finswap = finnow;
                finnow = finother; 
                finother = finswap;
            }

            int bytcalen = 0;
            double[] bytca = new double[8];
            if (bdytail != 0.0)
            {
                bytcalen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, ca, bdytail, bytca);
                var temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bytcalen, bytca, 2.0 * bdy,
                                                      temp16a);

                double[] bytcc = new double[8];
                var bytcclen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, cc, bdytail, bytcc);
                var temp16blen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bytcclen, bytcc, adx, temp16b);

                double[] bytaa = new double[8];
                var bytaalen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, aa, bdytail, bytaa);
                var temp16clen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bytaalen, bytaa, -cdx, temp16c);
                var temp32alen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                        temp16blen, temp16b, temp32a);
                var temp48len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16clen, temp16c,
                                                        temp32alen, temp32a, temp48);
                finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                        temp48, finother);
                finswap = finnow;
                finnow = finother; 
                finother = finswap;
            }

            int cxtablen = 0;
            double[] cxtab = new double[8];
            if (cdxtail != 0.0)
            {
                cxtablen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, ab, cdxtail, cxtab);
                var temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cxtablen, cxtab, 2.0 * cdx,
                                                      temp16a);
                double[] cxtbb = new double[8];
                var cxtbblen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, bb, cdxtail, cxtbb);
                var temp16blen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cxtbblen, cxtbb, ady, temp16b);
                double[] cxtaa = new double[8];
                var cxtaalen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, aa, cdxtail, cxtaa);
                var temp16clen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cxtaalen, cxtaa, -bdy, temp16c);

                var temp32alen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                        temp16blen, temp16b, temp32a);
                var temp48len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16clen, temp16c,
                                                        temp32alen, temp32a, temp48);
                finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                        temp48, finother);
                finswap = finnow;
                finnow = finother; 
                finother = finswap;
            }

            int cytablen = 0;
            double[] cytab = new double[8];
            if (cdytail != 0.0)
            {
                cytablen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, ab, cdytail, cytab);
                var temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cytablen, cytab, 2.0 * cdy,
                                                      temp16a);

                double[] cytaa = new double[8];
                var cytaalen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, aa, cdytail, cytaa);
                var temp16blen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cytaalen, cytaa, bdx, temp16b);

                double[] cytbb = new double[8];
                var cytbblen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, bb, cdytail, cytbb);
                var temp16clen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cytbblen, cytbb, -adx, temp16c);

                var temp32alen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                        temp16blen, temp16b, temp32a);
                var temp48len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16clen, temp16c,
                                                        temp32alen, temp32a, temp48);
                finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                        temp48, finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap;
            }

            int bcttlen;
            double[] bct = new double[8];
            int bctlen;
            var bctt = new double[4];
            double[] temp64 = new double[64];
            double[] u = new double[5];
            double[] v = new double[5];
            if ((adxtail != 0.0) || (adytail != 0.0))
            {
                if ((bdxtail != 0.0) || (bdytail != 0.0)
                    || (cdxtail != 0.0) || (cdytail != 0.0))
                {
                    MacrosHelpers.TwoProduct(bdxtail, cdy, out double ti1, out double ti0);
                    MacrosHelpers.TwoProduct(bdx, cdytail, out double tj1, out double tj0);
                    MacrosHelpers.TwoTwoSum(ti1, ti0, tj1, tj0, out u[3], out u[2], out u[1], out u[0]);
                    var negate = -bdy;
                    MacrosHelpers.TwoProduct(cdxtail, negate, out ti1, out ti0);
                    negate = -bdytail;
                    MacrosHelpers.TwoProduct(cdx, negate, out tj1, out tj0);

                    MacrosHelpers.TwoTwoSum(ti1, ti0, tj1, tj0, out v[3], out v[2], out v[1], out v[0]);
                    bctlen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(4, u, 4, v, bct);

                    MacrosHelpers.TwoProduct(bdxtail, cdytail, out ti1, out ti0);
                    MacrosHelpers.TwoProduct(cdxtail, bdytail, out tj1, out tj0);
                    MacrosHelpers.TwoTwoDiff(ti1, ti0, tj1, tj0, out bctt[3], out bctt[2], out bctt[1], out bctt[0]);
                    bcttlen = 4;
                }
                else
                {
                    bct[0] = 0.0;
                    bctlen = 1;
                    bctt[0] = 0.0;
                    bcttlen = 1;
                }

                double[] temp8 = new double[8];
                if (adxtail != 0.0)
                {
                    var temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(axtbclen, axtbc, adxtail, temp16a);
                    double[] axtbct = new double[16];
                    var axtbctlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bctlen, bct, adxtail, axtbct);
                    var temp32alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(axtbctlen, axtbct, 2.0 * adx,
                                                          temp32a);
                    var temp48len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp32alen, temp32a, temp48);
                    finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                            temp48, finother);
                    finswap = finnow; 
                    finnow = finother;
                    finother = finswap;
                    if (bdytail != 0.0)
                    {
                        var temp8len = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, cc, adxtail, temp8);
                        temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(temp8len, temp8, bdytail,
                                                              temp16a);
                        finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp16alen,
                                                                temp16a, finother);
                        finswap = finnow;
                        finnow = finother; 
                        finother = finswap;
                    }
                    if (cdytail != 0.0)
                    {
                        var temp8len = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, bb, -adxtail, temp8);
                        temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(temp8len, temp8, cdytail,
                                                              temp16a);
                        finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp16alen,
                                                                temp16a, finother);
                        finswap = finnow;
                        finnow = finother;
                        finother = finswap;
                    }

                    temp32alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(axtbctlen, axtbct, adxtail,
                                                          temp32a);
                    double[] axtbctt = new double[8];
                    var axtbcttlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bcttlen, bctt, adxtail, axtbctt);
                    temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(axtbcttlen, axtbctt, 2.0 * adx,
                                                          temp16a);
                    var temp16blen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(axtbcttlen, axtbctt, adxtail,
                                                          temp16b);
                    var temp32blen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp16blen, temp16b, temp32b);
                    var temp64len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp32alen, temp32a,
                                                            temp32blen, temp32b, temp64);
                    finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp64len,
                                                            temp64, finother);
                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;
                }
                if (adytail != 0.0)
                {
                    double[] aytbct = new double[16];
                    var temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(aytbclen, aytbc, adytail, temp16a);
                    var aytbctlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bctlen, bct, adytail, aytbct);
                    var temp32alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(aytbctlen, aytbct, 2.0 * ady,
                                                          temp32a);
                    var temp48len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp32alen, temp32a, temp48);
                    finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                            temp48, finother);

                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;

                    temp32alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(aytbctlen, aytbct, adytail,
                                                          temp32a);
                    double[] aytbctt = new double[16];
                    var aytbcttlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bcttlen, bctt, adytail, aytbctt);
                    temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(aytbcttlen, aytbctt, 2.0 * ady,
                                                          temp16a);
                    var temp16blen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(aytbcttlen, aytbctt, adytail,
                                                          temp16b);
                    var temp32blen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp16blen, temp16b, temp32b);
                    var temp64len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp32alen, temp32a,
                                                            temp32blen, temp32b, temp64);
                    finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp64len,
                                                            temp64, finother);
                    finswap = finnow; 
                    finnow = finother;
                    finother = finswap;
                }
            }
            double[] cat = new double[8];
            double[] catt = new double[4];
            if ((bdxtail != 0.0) || (bdytail != 0.0))
            {
                int cattlen;
                int catlen;
                if ((cdxtail != 0.0) || (cdytail != 0.0)
                    || (adxtail != 0.0) || (adytail != 0.0))
                {
                    MacrosHelpers.TwoProduct(cdxtail, ady, out double ti1, out double ti0);
                    MacrosHelpers.TwoProduct(cdx, adytail, out double tj1, out double tj0);

                    MacrosHelpers.TwoTwoSum(ti1, ti0, tj1, tj0, out u[3], out u[2], out u[1], out u[0]);

                    var negate = -cdy;
                    MacrosHelpers.TwoProduct(adxtail, negate, out ti1, out ti0);
                    negate = -cdytail;
                    MacrosHelpers.TwoProduct(adx, negate, out tj1, out tj0);
                    MacrosHelpers.TwoTwoSum(ti1, ti0, tj1, tj0, out v[3], out v[2], out v[1], out v[0]);

                    catlen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(4, u, 4, v, cat);

                    MacrosHelpers.TwoProduct(cdxtail, adytail, out ti1, out ti0);
                    MacrosHelpers.TwoProduct(adxtail, cdytail, out tj1, out tj0);
                    MacrosHelpers.TwoTwoDiff(ti1, ti0, tj1, tj0, out catt[3], out catt[2], out catt[1], out catt[0]);
                    cattlen = 4;
                }
                else
                {
                    cat[0] = 0.0;
                    catlen = 1;
                    catt[0] = 0.0;
                    cattlen = 1;
                }

                double[] temp8 = new double[8];
                double[] bxtcat = new double[16];
                double[] bxtcatt = new double[8];
                if (bdxtail != 0.0)
                {
                    var temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bxtcalen, bxtca, bdxtail, temp16a);
                    var bxtcatlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(catlen, cat, bdxtail, bxtcat);
                    var temp32alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bxtcatlen, bxtcat, 2.0 * bdx,
                                                          temp32a);
                    var temp48len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                             temp32alen, temp32a, temp48);
                    finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                            temp48, finother);
                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;
                    if (cdytail != 0.0)
                    {
                        var temp8len = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, aa, bdxtail, temp8);
                        temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(temp8len, temp8, cdytail,
                                                              temp16a);
                        finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp16alen,
                                                                temp16a, finother);
                        finswap = finnow;
                        finnow = finother;
                        finother = finswap;
                    }
                    if (adytail != 0.0)
                    {
                        var temp8len = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, cc, -bdxtail, temp8);
                        temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(temp8len, temp8, adytail,
                                                              temp16a);
                        finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp16alen,
                                                                temp16a, finother);
                        finswap = finnow; 
                        finnow = finother;
                        finother = finswap;
                    }

                    temp32alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bxtcatlen, bxtcat, bdxtail,
                                                          temp32a);
                    var bxtcattlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cattlen, catt, bdxtail, bxtcatt);
                    temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bxtcattlen, bxtcatt, 2.0 * bdx,
                                                          temp16a);
                    var temp16blen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bxtcattlen, bxtcatt, bdxtail,
                                                          temp16b);
                    var temp32blen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp16blen, temp16b, temp32b);
                    var temp64len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp32alen, temp32a,
                                                            temp32blen, temp32b, temp64);
                    finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp64len,
                                                            temp64, finother);
                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;
                }
                double[] bytcat = new double[16];
                double[] bytcatt = new double[8];
                if (bdytail != 0.0)
                {
                    var temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bytcalen, bytca, bdytail, temp16a);
                    var bytcatlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(catlen, cat, bdytail, bytcat);
                    var temp32alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bytcatlen, bytcat, 2.0 * bdy,
                                                          temp32a);
                    var temp48len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp32alen, temp32a, temp48);
                    finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                            temp48, finother);
                    finswap = finnow;
                    finnow = finother; 
                    finother = finswap;


                    temp32alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bytcatlen, bytcat, bdytail,
                                                          temp32a);
                    var bytcattlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cattlen, catt, bdytail, bytcatt);
                    temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bytcattlen, bytcatt, 2.0 * bdy,
                                                          temp16a);
                    var temp16blen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(bytcattlen, bytcatt, bdytail,
                                                          temp16b);
                    var temp32blen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp16blen, temp16b, temp32b);
                    var temp64len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp32alen, temp32a,
                                                            temp32blen, temp32b, temp64);
                    finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp64len,
                                                            temp64, finother);
                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;
                }
            }
            double[] abt = new double[8];
            double[] abtt = new double[4];
            if ((cdxtail != 0.0) || (cdytail != 0.0))
            {
                int abtlen;
                int abttlen;
                if ((adxtail != 0.0) || (adytail != 0.0)
                    || (bdxtail != 0.0) || (bdytail != 0.0))
                {
                    MacrosHelpers.TwoProduct(adxtail, bdy, out double ti1, out double ti0);
                    MacrosHelpers.TwoProduct(adx, bdytail, out double tj1, out double tj0);
                    MacrosHelpers.TwoTwoSum(ti1, ti0, tj1, tj0, out u[3], out u[2], out u[1], out u[0]);
                    var negate = -ady;
                    MacrosHelpers.TwoProduct(bdxtail, negate, out ti1, out ti0);
                    negate = -adytail;
                    MacrosHelpers.TwoProduct(bdx, negate, out tj1, out tj0);
                    MacrosHelpers.TwoTwoSum(ti1, ti0, tj1, tj0, out v[3], out v[2], out v[1], out v[0]);
                    abtlen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(4, u, 4, v, abt);

                    MacrosHelpers.TwoProduct(adxtail, bdytail, out ti1, out ti0);
                    MacrosHelpers.TwoProduct(bdxtail, adytail, out tj1, out tj0);
                    MacrosHelpers.TwoTwoDiff(ti1, ti0, tj1, tj0, out abtt[3], out abtt[2], out abtt[1], out abtt[0]);
                    abttlen = 4;
                }
                else
                {
                    abt[0] = 0.0;
                    abtlen = 1;
                    abtt[0] = 0.0;
                    abttlen = 1;
                }

                double[] temp8 = new double[8];
                double[] cxtabt = new double[16];
                double[] cxtabtt = new double[8];
                if (cdxtail != 0.0)
                {
                    var temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cxtablen, cxtab, cdxtail, temp16a);
                    var cxtabtlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(abtlen, abt, cdxtail, cxtabt);
                    var temp32alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cxtabtlen, cxtabt, 2.0 * cdx,
                                                          temp32a);
                    var temp48len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp32alen, temp32a, temp48);
                    finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                            temp48, finother);
                    finswap = finnow;
                    finnow = finother;
                    finother = finswap;

                    if (adytail != 0.0)
                    {
                        var temp8len = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, bb, cdxtail, temp8);
                        temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(temp8len, temp8, adytail,
                                                              temp16a);
                        finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp16alen,
                                                                temp16a, finother);
                        finswap = finnow; 
                        finnow = finother; 
                        finother = finswap;
                    }
                    if (bdytail != 0.0)
                    {
                        var temp8len = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(4, aa, -cdxtail, temp8);
                        temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(temp8len, temp8, bdytail,
                                                              temp16a);
                        finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp16alen,
                                                                temp16a, finother);
                        finswap = finnow; 
                        finnow = finother; 
                        finother = finswap;
                    }

                    temp32alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cxtabtlen, cxtabt, cdxtail,
                                                          temp32a);
                    var cxtabttlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(abttlen, abtt, cdxtail, cxtabtt);
                    temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cxtabttlen, cxtabtt, 2.0 * cdx,
                                                          temp16a);
                    var temp16blen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cxtabttlen, cxtabtt, cdxtail,
                                                          temp16b);
                    var temp32blen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp16blen, temp16b, temp32b);
                    var temp64len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp32alen, temp32a,
                                                            temp32blen, temp32b, temp64);
                    finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp64len,
                                                            temp64, finother);
                    finswap = finnow; 
                    finnow = finother; 
                    finother = finswap;
                }

                double[] cytabt = new double[16];
                double[] cytabtt = new double[8];
                if (cdytail != 0.0)
                {
                    var temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cytablen, cytab, cdytail, temp16a);
                    var cytabtlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(abtlen, abt, cdytail, cytabt);
                    var temp32alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cytabtlen, cytabt, 2.0 * cdy,
                                                          temp32a);
                    var temp48len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp32alen, temp32a, temp48);
                    finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp48len,
                                                            temp48, finother);
                    finswap = finnow; finnow = finother; finother = finswap;


                    temp32alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cytabtlen, cytabt, cdytail,
                                                          temp32a);
                    var cytabttlen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(abttlen, abtt, cdytail, cytabtt);
                    temp16alen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cytabttlen, cytabtt, 2.0 * cdy,
                                                          temp16a);
                    var temp16blen = ArithmeticFunctionsHelpers.ScaleExpansionZeroelim(cytabttlen, cytabtt, cdytail,
                                                          temp16b);
                    var temp32blen = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp16alen, temp16a,
                                                            temp16blen, temp16b, temp32b);
                    var temp64len = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(temp32alen, temp32a,
                                                            temp32blen, temp32b, temp64);
                    finlength = ArithmeticFunctionsHelpers.FastExpansionSumZeroelim(finlength, finnow, temp64len,
                                                            temp64, finother);
                    finswap = finnow; 
                    finnow = finother; 
                    finother = finswap;
                }
            }

            return finnow[finlength - 1];
        }

        public static double Robust(double[] pa, double[] pb, double[] pc, double[] pd)
        {
            var adx = pa[0] - pd[0];
            var bdx = pb[0] - pd[0];
            var cdx = pc[0] - pd[0];
            var ady = pa[1] - pd[1];
            var bdy = pb[1] - pd[1];
            var cdy = pc[1] - pd[1];

            var bdxcdy = bdx * cdy;
            var cdxbdy = cdx * bdy;
            var alift = adx * adx + ady * ady;

            var cdxady = cdx * ady;
            var adxcdy = adx * cdy;
            var blift = bdx * bdx + bdy * bdy;

            var adxbdy = adx * bdy;
            var bdxady = bdx * ady;
            var clift = cdx * cdx + cdy * cdy;

            var det = alift * (bdxcdy - cdxbdy)
                  + blift * (cdxady - adxcdy)
                  + clift * (adxbdy - bdxady);

            var permanent = (Math.Abs(bdxcdy) + Math.Abs(cdxbdy)) * alift
                        + (Math.Abs(cdxady) + Math.Abs(adxcdy)) * blift
                        + (Math.Abs(adxbdy) + Math.Abs(bdxady)) * clift;

            var errbound = MacrosHelpers.IccerrboundA * permanent;
            if ((det > errbound) || (-det > errbound))
            {
                return det;
            }

            return Adapt(pa, pb, pc, pd, permanent);
        }
    }
}
