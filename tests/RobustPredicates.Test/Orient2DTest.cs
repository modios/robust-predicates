using System;
using System.IO;
using System.Linq;
using Xunit;

namespace RobustPredicates.Test
{
    public class Orient2DTest
    {
        private const int NSimpleData = 100;
        private const int NPscicoData = 1000;

        [Fact]
        public void Robust_ShouldSucceed()
        {
            Assert.True(Orient2D.Robust
                (new double[] { 0, 0 }, 
                new double[] { 1, 1 },
                new double[] { 0, 1 }) > 0);
            Assert.True(Orient2D.Robust(
                new double[] { 0, 0 }, 
                new double[] { 0, 1 },
                new double[] { 1, 1 }) < 0);
            Assert.True(Orient2D.Robust(
                new double[] { 0, 0 },
                new double[] { 0.5, 0.5 },
                new double[] { 1, 1 }) == 0);
        }

        [Fact]
        public void Robust_FromFile_ShouldSucceed()
        {
            double[] points =
              File.ReadAllLines("test_data/simple_data/orient2d-points2d.txt")
              .Select(n => n.Split()).SelectMany(x => x)
              .Select(s => double.Parse(s)).ToArray();
            double[] results =
             File.ReadAllLines("test_data/simple_data/results-orient2d.txt")
             .Select(n => n.Split()).SelectMany(x => x)
             .Select(s => double.Parse(s)).ToArray();

            int count = 0;
            for (int i = 0; i < NSimpleData; i += 6)
            {
                Assert.Equal(Math.Sign(results[count++]), Math.Sign(Orient2D.Robust(
                    new double[] { points[i], points[i + 1] },
                    new double[] { points[i + 2], points[i + 3] },
                    new double[] { points[i + 4], points[i + 5] })));
            }
        }

        [Fact]
        public void Fast_ShouldSucceed()
        {
            Assert.True(Orient2D.Fast(
                new double[] { 0, 0 }, 
                new double[] { 1, 1 },
                new double[] { 0, 1 }) > 0);
            Assert.True(Orient2D.Fast(
                new double[] { 0, 0 }, 
                new double[] { 0, 1 },
                new double[] { 1, 1 }) < 0);
            Assert.True(Orient2D.Fast(
                new double[] { 0, 0 },
                new double[] { 0.5, 0.5 },
                new double[] { 1, 1 }) == 0);
        }

        [Fact]
        public void Fast_FromFile_ShouldSucceed()
        {
            double[] points =
              File.ReadAllLines("test_data/simple_data/orient2d-points2d.txt")
              .Select(n => n.Split()).SelectMany(x => x)
              .Select(s => double.Parse(s)).ToArray();
            double[] results =
             File.ReadAllLines("test_data/simple_data/results-orient2d.txt")
             .Select(n => n.Split()).SelectMany(x => x)
             .Select(s => double.Parse(s)).ToArray();

            int count = 0;
            for (int i = 0; i < NSimpleData; i += 6)
            {
                Assert.Equal(results[count++], Orient2D.Fast(
                    new double[] { points[i], points[i + 1] },
                    new double[] { points[i + 2], points[i + 3] },
                    new double[] { points[i + 4], points[i + 5] }), 4);
            }
        }

        [Fact]
        public void Robust_FromFile_Pscico_ShouldSucceed()
        {
            double[] points =
              File.ReadAllLines("test_data/pscico_data/orient2d.txt")
              .Select(n => n.Split()).SelectMany(x => x)
              .Select(s => double.Parse(s)).ToArray();

            for (int i = 0; i < NPscicoData; i += 8)
            {
                Assert.Equal(Math.Sign(points[i + 7]), Math.Sign(Orient2D.Robust(
                    new double[] { points[i + 1], points[i + 2] },
                    new double[] { points[i + 3], points[i + 4] },
                    new double[] { points[i + 5], points[i + 6] })));
            }
        }
    }
}
