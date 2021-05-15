using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RobustPredicates.Test
{
    public static class TestHelpers
    {
        private static Random _random = default(Random);

        private static double GetRandomDouble(double min, double max) =>   (max - min) * _random.NextDouble() + min;
                
        public static double[] Create2DRandomPoints(double min, double max, int seed = 0)
        {
            if(_random == default(Random))
            {
                _random = new Random(seed);
            }

            GetRandomDouble(min, max);

            return Enumerable.Range(0, 6).Select(s => GetRandomDouble(min, max)).ToArray();
        }

        public static async void Create2DPointsAndWriteToFile(int numberOfPoint, double min, double max)
        {
            var lines = Enumerable.Range(0, numberOfPoint).Select(p => Create2DRandomPoints(min, max)).Select(
                r => $"{r[0]} {r[1]} {r[2]} {r[3]} {r[4]} {r[5]}" );
            await File.WriteAllLinesAsync("points2D.txt", lines);
        }
    }
}
