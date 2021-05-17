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
                
        public static double[] Create2DRandomPoints(double min, double max, int numberOfPoints = 3, int seed = 0)
        {
            if(_random == default(Random))
            {
                _random = new Random(seed);
            }

            GetRandomDouble(min, max);

            return Enumerable.Range(0, numberOfPoints * 2).Select(s => GetRandomDouble(min, max)).ToArray();
        }

        public static async void Create2DPointsAndWriteToFile(int numberOfPoint, double min, double max, int perRow = 3 , string filename ="points2D.txt")
        {
            var lines = Enumerable.Range(0, numberOfPoint).Select(p => Create2DRandomPoints(min, max, perRow)).Select(
                r => {
                    var stringBuilder = new StringBuilder();
                    stringBuilder.Append(r[0]);
                    for(int i=1; i< r.Length; i++)
                    {
                        stringBuilder.Append(" " + r[i] );
                    }

                    return stringBuilder.ToString();
                }).ToArray();

            await File.WriteAllLinesAsync(filename, lines);
        }
    }
}
