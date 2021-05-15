namespace RobustPredicates
{
    internal static class MacrosHelpers
    {
        public static double Epsilon;
        public static double Splitter;
        public static double Resulterrbound;
        public static double CcwerrboundA;
        public static double CcwerrboundB;
        public static double CcwerrboundC;
        public static double O3derrboundA;
        public static double O3derrboundB;
        public static double O3derrboundC;
        public static double IccerrboundA;
        public static double IccerrboundB;
        public static double IccerrboundC;
        public static double IsperrboundA;
        public static double IsperrboundB;
        public static double IsperrboundC;

        static MacrosHelpers()
        {
            double half;
            double check, lastcheck;

            bool everyOther = true;

            half = 0.5;
            Epsilon = 1.0;
            Splitter = 1.0;
            check = 1.0;

            do
            {
                lastcheck = check;
                Epsilon *= half;
                if (everyOther)
                {
                    Splitter *= 2.0;
                }
                everyOther = !everyOther;
                check = 1.0 + Epsilon;
            } while ((check != 1.0) && (check != lastcheck));

            Splitter += 1.0;

            Resulterrbound = (3.0 + 8.0 * Epsilon) * Epsilon;
            CcwerrboundA = (3.0 + 16.0 * Epsilon) * Epsilon;
            CcwerrboundB = (2.0 + 12.0 * Epsilon) * Epsilon;
            CcwerrboundC = (9.0 + 64.0 * Epsilon) * Epsilon * Epsilon;
            O3derrboundA = (7.0 + 56.0 * Epsilon) * Epsilon;
            O3derrboundB = (3.0 + 28.0 * Epsilon) * Epsilon;
            O3derrboundC = (26.0 + 288.0 * Epsilon) * Epsilon * Epsilon;
            IccerrboundA = (10.0 + 96.0 * Epsilon) * Epsilon;
            IccerrboundB = (4.0 + 48.0 * Epsilon) * Epsilon;
            IccerrboundC = (44.0 + 576.0 * Epsilon) * Epsilon * Epsilon;
            IsperrboundA = (16.0 + 224.0 * Epsilon) * Epsilon;
            IsperrboundB = (5.0 + 72.0 * Epsilon) * Epsilon;
            IsperrboundC = (71.0 + 1408.0 * Epsilon) * Epsilon * Epsilon;
        }

        internal static void FastTwoSumTail(double a, double b, double x, out double y)
        {
            double bvirt = x - a;
            y = b - bvirt;
        }

        internal static void FastTwoSum(double a, double b, out double x, out double y)
        {
            x = a + b;
            FastTwoSumTail(a, b, x, out y);
        }

        internal static void FastTwoDiffTail(double a, double b, double x, out double y)
        {
            double bvirt = a - x;
            y = bvirt - b;
        }

        internal static void FastTwoDiff(double a, double b, out double x, out double y)
        {
            x = a - b;
            FastTwoDiffTail(a, b, x, out y);
        }

        internal static void TwoSumTail(double a, double b, double x, out double y)
        {
            double bvirt = x - a;
            double avirt = x - bvirt;
            double bround = b - bvirt;
            double around = a - avirt;
            y = around + bround;
        }

        internal static void TwoSum(double a, double b, out double x, out double y)
        {
            x = a + b;
            TwoSumTail(a, b, x, out y);
        }

        internal static void TwoDiffTail(double a, double b, double x, out double y)
        {
            double bvirt = a - x;
            double avirt = x + bvirt;
            double bround = bvirt - b;
            double around = a - avirt;
            y = around + bround;
        }

        internal static void TwoDiff(double a, double b, out double x, out double y)
        {
            x = a - b;
            TwoDiffTail(a, b, x, out y);
        }

        internal static void Split(double a, out double ahi, out double alo)
        {
            double c = (double)(Splitter * a);
            double abig = (double)(c - a);
            ahi = c - abig;
            alo = a - ahi;
        }

        internal static void TwoProductTail(double a, double b, double x, out double y)
        {
            Split(a, out double ahi, out double alo);
            Split(b, out double bhi, out double blo);
            double err1 = x - (ahi * bhi);
            double err2 = err1 - (alo * bhi);
            double err3 = err2 - (ahi * blo);
            y = (alo * blo) - err3;
        }

        internal static void TwoProduct(double a, double b, out double x, out double y)
        {
            x = a * b;
            TwoProductTail(a, b, x, out y);
        }

        internal static void TwoProductPresplit(double a, double b, double bhi, double blo, out double x, out double y)
        {
            x = a * b;
            Split(a, out double ahi, out double alo);
            double err1 = x - (ahi * bhi);
            double err2 = err1 - (alo * bhi);
            double err3 = err2 - (ahi * blo);
            y = (alo * blo) - err3;
        }

        internal static void TwoProduct2Presplit(double a, double ahi, double alo, double b, double bhi, double blo, out double x, out double y)
        {
            x = a * b;
            double err1 = x - (ahi * bhi);
            double err2 = err1 - (alo * bhi);
            double err3 = err2 - (ahi * blo);
            y = (alo * blo) - err3;
        }

        internal static void SquareTail(double a, double x, out double y)
        {
            Split(a, out double ahi, out double alo);
            double err1 = x - (ahi * ahi);
            double err3 = err1 - ((ahi + ahi) * alo);
            y = (alo * alo) - err3;
        }

        internal static void Square(double a, out double x, out double y)
        {
            x = a * a;
            SquareTail(a, x, out y);
        }

        internal static void TwoOneSum(double a1, double a0, double b, out double x2, out double x1, out double x0)
        {
            TwoSum(a0, b, out double _i, out x0);
            TwoSum(a1, _i, out x2, out x1);
        }

        internal static void TwoOneDiff(double a1, double a0, double b, out double x2, out double x1, out double x0)
        {
            TwoDiff(a0, b, out double _i, out x0);
            TwoSum(a1, _i, out x2, out x1);
        }

        internal static void TwoTwoSum(double a1, double a0, double b1, double b0, out double x3, out double x2, out double x1, out double x0)
        {
            TwoOneSum(a1, a0, b0, out double _j, out double _0, out x0);
            TwoOneSum(_j, _0, b1, out x3, out x2, out x1);
        }

        internal static void TwoTwoDiff(double a1, double a0, double b1, double b0, out double x3, out double x2, out double x1, out double x0)
        {
            TwoOneDiff(a1, a0, b0, out double _j, out double _0, out x0);
            TwoOneDiff(_j, _0, b1, out x3, out x2, out x1);
        }

        internal static void FourOneSum(double a3, double a2, double a1, double a0, double b, out double x4, out double x3, out double x2, out double x1, out double x0)
        {
            TwoOneSum(a1, a0, b, out double _j, out x1, out x0);
            TwoOneSum(a3, a2, _j, out x4, out x3, out x2);
        }

        internal static void FourTwoSum(double a3, double a2, double a1, double a0, double b1, double b0, out double x5, out double x4, out double x3, out double x2, out double x1, out double x0)
        {
            FourOneSum(a3, a2, a1, a0, b0, out double _k, out double _2, out double _1, out double _0, out x0);
            FourOneSum(_k, _2, _1, _0, b1, out x5, out x4, out x3, out x2, out x1);
        }

        internal static void FourFourSum(double a3, double a2, double a1, double a0, double b4, double b3, double b1, double b0,
                                          out double x7, out double x6, out double x5, out double x4, out double x3, out double x2,
                                          out double x1, out double x0)
        {
            FourTwoSum(a3, a2, a1, a0, b1, b0, out double _l, out double _2, out double _1, out double _0, out x1, out x0);
            FourTwoSum(_l, _2, _1, _0, b4, b3, out x7, out x6, out x5, out x4, out x3, out x2);
        }

        internal static void EightOneSum(double a7, double a6, double a5, double a4, double a3, double a2, double a1, double a0, double b,
            out double x8, out double x7, out double x6, out double x5, out double x4, out double x3, out double x2, out double x1, out double x0)
        {
            FourOneSum(a3, a2, a1, a0, b, out double _j, out x3, out x2, out x1, out x0);
            FourOneSum(a7, a6, a5, a4, _j, out x8, out x7, out x6, out x5, out x4);
        }

        internal static void EightTwoSum(double a7, double a6, double a5, double a4, double a3, double a2, double a1, double a0, double b1, double b0,
            out double x9, out double x8, out double x7, out double x6, out double x5, out double x4, out double x3, out double x2, out double x1, out double x0)
        {
            EightOneSum(a7, a6, a5, a4, a3, a2, a1, a0, b0, out double _k, out double _6, out double _5, out double _4, out double _3, out double _2, out double _1, out double _0, out x0);
            EightOneSum(_k, _6, _5, _4, _3, _2, _1, _0, b1, out x9, out x8, out x7, out x6, out x5, out x4, out x3, out x2, out x1);
        }

        internal static void TwoTwoProduct(double a1, double a0, double b1, double b0,
            out double x7, out double x6, out double x5, out double x4, out double x3, out double x2, out double x1, out double x0)
        {
            Split(a0, out double a0hi, out double a0lo);
            Split(b0, out double bhi, out double blo);
            TwoProduct2Presplit(a0, a0hi, a0lo, b0, bhi, blo, out double _i, out x0);
            Split(a1, out double a1hi, out double a1lo);
            TwoProduct2Presplit(a1, a1hi, a1lo, b0, bhi, blo, out double _j, out double _0);
            TwoSum(_i, _0, out double _k, out double _1);
            FastTwoSum(_j, _k, out double _l, out double _2);
            Split(b1, out bhi, out blo);
            TwoProduct2Presplit(a0, a0hi, a0lo, b1, bhi, blo, out _i, out _0);
            TwoSum(_1, _0, out _k, out x1);
            TwoSum(_2, _k, out _j, out _1);
            TwoSum(_l, _j, out double _m, out _2);
            TwoProduct2Presplit(a1, a1hi, a1lo, b1, bhi, blo, out _j, out _0);
            TwoSum(_i, _0, out double _n, out _0);
            TwoSum(_1, _0, out _i, out x2);
            TwoSum(_2, _i, out _k, out _1);
            TwoSum(_m, _k, out _l, out _2);
            TwoSum(_j, _n, out _k, out _0);
            TwoSum(_1, _0, out _j, out x3);
            TwoSum(_2, _j, out _i, out _1);
            TwoSum(_l, _i, out _m, out _2);
            TwoSum(_1, _k, out _i, out x4);
            TwoSum(_2, _i, out _k, out x5);
            TwoSum(_m, _k, out x7, out x6);
        }
    }
}
