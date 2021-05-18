using System;
using System.IO;
using System.Linq;
using Xunit;

namespace RobustPredicates.Test
{
    public class Orient3DTest
    {
        [Fact]
        public void Robust_ShouldSucceed()
        {
            Assert.True(Orient3D.Robust(new double[] { 0, 0, 0 }, new double[] { 1, 0, 0 }, new double[] { 0 , 1, 1 }, new double[] { 0, 0, 1 }) < 0);
            Assert.True(Orient3D.Robust(new double[] { 0, 0, 0 }, new double[] { 0, 1, 0 }, new double[] { 1, 0, 0 }, new double[] { 0, 0, 1 }) > 0);
            Assert.True(Orient3D.Robust(new double[] { 0, 0, 0 }, new double[] { 1, 1, 0 }, new double[] { 0, 1, 1 }, new double[] { 0, 1, 1 }) == 0);
        }

        [Fact]
        public void Robust_FromFile_ShouldSucceed()
        {
            double[] numbers =
              File.ReadAllLines("testData/oient3Dpoints.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();
            double[] results =
             File.ReadAllLines("testData/orient3DResults.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();
            int count = 0;
            for (int i = 0; i < 1000; i += 12)
            { 
                var orient3DResult = Orient3D.Robust(
                  new double[] { numbers[i], numbers[i + 1], numbers[i + 2] },
                  new double[] { numbers[i + 3], numbers[i + 4], numbers[i + 5] },
                  new double[] { numbers[i + 6], numbers[i + 7], numbers[i + 8] },
                  new double[] { numbers[i + 9], numbers[i + 10], numbers[i + 11] });
                Assert.Equal(Math.Sign(results[count++]), Math.Sign(orient3DResult));
            }
        }

        [Fact]
        public void Fast_ShouldSucceed()
        {
            Assert.True(Orient3D.Fast(new double[] { 0, 0, 0 }, new double[] { 1, 0, 0 }, new double[] { 0, 1, 1 }, new double[] { 0, 0, 1 }) < 0);
            Assert.True(Orient3D.Fast(new double[] { 0, 0, 0 }, new double[] { 0, 1, 0 }, new double[] { 1, 0, 0 }, new double[] { 0, 0, 1 }) > 0);
            Assert.True(Orient3D.Fast(new double[] { 0, 0, 0 }, new double[] { 1, 1, 0 }, new double[] { 0, 1, 1 }, new double[] { 0, 1, 1 }) == 0);
        }

        [Fact]
        public void Fast_FromFile_ShouldSucceed()
        {
            double[] numbers =
              File.ReadAllLines("testData/oient3Dpoints.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();
            double[] results =
             File.ReadAllLines("testData/orient3DResults.txt").Select(n => n.Split()).SelectMany(x => x).Select(s => double.Parse(s)).ToArray();
            int count = 0;
            for (int i = 0; i < 1000; i += 12)
            {
                var orient3DResult = Orient3D.Fast(
                  new double[] { numbers[i], numbers[i + 1], numbers[i + 2] },
                  new double[] { numbers[i + 3], numbers[i + 4], numbers[i + 5] },
                  new double[] { numbers[i + 6], numbers[i + 7], numbers[i + 8] },
                  new double[] { numbers[i + 9], numbers[i + 10], numbers[i + 11] });
                Assert.Equal(Math.Sign(results[count++]), Math.Sign(orient3DResult));
            }
        }
    }
}
