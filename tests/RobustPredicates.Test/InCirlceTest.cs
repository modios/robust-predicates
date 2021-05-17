using System.IO;
using System.Linq;
using Xunit;

namespace RobustPredicates.Test
{
    public class InCirlceTest
    {
        [Fact]
        public void Fast_ShouldSucceed()
        {
            Assert.True(InCirlce.Fast(new double[] { 0, -1 }, new double[] { 0, 1 }, new double[] { 1, 0 }, new double[] { 0, 0 }) < 0);
            Assert.True(InCirlce.Fast(new double[] { 0, -1 }, new double[] { 0, 1 }, new double[] { 1, 0 }, new double[] { -1, 0 }) == 0);
            Assert.True(InCirlce.Fast(new double[] { 0, -1 }, new double[] { 0, 1}, new double[] { 1, 0, }, new double[] { 10, 10 }) > 0);
        }

        [Fact]
        public void Fast_FromFile_ShouldSucceed()
        {
            double[] numbers =
              File.ReadAllLines("testData/inCriclePoinst2D.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();
            double[] results =
             File.ReadAllLines("testData/inCircleResults.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();
            int count = 0;
            for (int i = 0; i < 100; i += 8)
            {
                Assert.Equal(results[count++], InCirlce.Fast(
                    new double[] { numbers[i], numbers[i + 1] },
                    new double[] { numbers[i + 2], numbers[i + 3] },
                    new double[] { numbers[i + 4], numbers[i + 5] },
                    new double[] { numbers[i + 6], numbers[i + 7] }), 5);
            }
        }

        [Fact]
        public void Robust_ShouldSucceed()
        {
            Assert.True(InCirlce.Robust(new double[] { 0, -1 }, new double[] { 0, 1 }, new double[] { 1, 0 }, new double[] { 0, 0 }) < 0);
            Assert.True(InCirlce.Robust(new double[] { 0, -1 }, new double[] { 0, 1 }, new double[] { 1, 0 }, new double[] { -1, 0 }) == 0);
            Assert.True(InCirlce.Robust(new double[] { 0, -1 }, new double[] { 0, 1 }, new double[] { 1, 0, }, new double[] { 10, 10 }) > 0);
        }
        
        [Fact]
        public void Robust_FromFile_Set1_ShouldSucceed()
        {
            double[] numbers =
              File.ReadAllLines("testData/inCriclePoinst2D.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();
            double[] results =
             File.ReadAllLines("testData/inCircleResults.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();
            int count = 0;
            for (int i = 0; i < 100; i += 8)
            {
                Assert.Equal(results[count++], InCirlce.Robust(
                    new double[] { numbers[i], numbers[i + 1] },
                    new double[] { numbers[i + 2], numbers[i + 3] },
                    new double[] { numbers[i + 4], numbers[i + 5] },
                    new double[] { numbers[i + 6], numbers[i + 7] }), 5);
            }
        }

        [Fact]
        public void Robust_FromFile_Set2_ShouldSucceed()
        {
            double[] numbers =
              File.ReadAllLines("testData/inCircle2d.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();

            for (int i = 0; i < 1000; i += 10)
            {
                Assert.Equal(numbers[i + 9], InCirlce.Robust(
                    new double[] { numbers[i + 1], numbers[i + 2] },
                    new double[] { numbers[i + 3], numbers[i + 4] },
                    new double[] { numbers[i + 5], numbers[i + 6] },
                    new double[] { numbers[i + 7], numbers[i + 8] }) > 0 ? 1 : -1);
            }
        }
    }
}
