using System.IO;
using System.Linq;
using Xunit;

namespace RobustPredicates.Test
{
    public class Orient2DTest
    {
        [Fact]
        public void Robust_ShouldSucceed()
        {
            Assert.True(Orient2D.Robust(new double[] { 0, 0 }, new double[] { 1, 1 }, new double[] { 0, 1 }) > 0);
            Assert.True(Orient2D.Robust(new double[] { 0, 0 }, new double[] { 0, 1 }, new double[] { 1, 1 }) < 0);
            Assert.True(Orient2D.Robust(new double[] { 0, 0 }, new double[] { 0.5, 0.5 }, new double[] { 1, 1 }) == 0);
        }

        [Fact]
        public void Robust_FromFile_ShouldSucceed()
        {
            double[] numbers =
              File.ReadAllLines("testData/points2D.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();
            double[] results =
             File.ReadAllLines("testData/orient2DResults.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();
            int count = 0;
            for (int i = 0; i < 100; i += 6)
            {
                Assert.Equal(results[count++], Orient2D.Robust(
                    new double[] { numbers[i], numbers[i + 1] },
                    new double[] { numbers[i + 2], numbers[i + 3] },
                    new double[] { numbers[i + 4], numbers[i + 5] }), 4);
            }
        }

        [Fact]
        public void Fast_ShouldSucceed()
        {
            Assert.True(Orient2D.Fast(new double[] { 0, 0 }, new double[] { 1, 1 }, new double[] { 0, 1 }) > 0);
            Assert.True(Orient2D.Fast(new double[] { 0, 0 }, new double[] { 0, 1 }, new double[] { 1, 1 }) < 0);
            Assert.True(Orient2D.Fast(new double[] { 0, 0 }, new double[] { 0.5, 0.5 }, new double[] { 1, 1 }) == 0);
        }

        [Fact]
        public void Fast_FromFile_ShouldSucceed()
        {
            double[] numbers =
              File.ReadAllLines("testData/points2D.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();
            double[] results =
             File.ReadAllLines("testData/orient2DResults.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();
            int count = 0;
            for (int i = 0; i < 100; i += 6)
            {
                Assert.Equal(results[count++], Orient2D.Fast(
                    new double[] { numbers[i], numbers[i + 1] },
                    new double[] { numbers[i + 2], numbers[i + 3] },
                    new double[] { numbers[i + 4], numbers[i + 5] }), 4);
            }
        }
    }
}
