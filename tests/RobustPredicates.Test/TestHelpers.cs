using System;
using System.IO;
using System.Linq;
using System.Text;

namespace RobustPredicates.Test
{
    public static class TestHelpers
    {
        private static Random _random = default(Random);

        private static double GetRandomDouble(double min, double max) => (max - min) * _random.NextDouble() + min;

        public static double[] CreateRandomPoints(double min, double max, int numberOfPoints = 3, int dimensions = 2, int seed = 0)
        {
            if (_random == default(Random))
            {
                _random = new Random(seed);
            }

            GetRandomDouble(min, max);

            return Enumerable.Range(0, numberOfPoints * dimensions).Select(s => GetRandomDouble(min, max)).ToArray();
        }

        public static async void CreatePointsAndWriteToFile(int numberOfPoint, double min, double max, int perRow = 3, int dimensions = 2, string filename = "points2D.txt")
        {
            var lines = Enumerable.Range(0, numberOfPoint).Select(p => CreateRandomPoints(min, max, perRow, dimensions)).Select(
                r =>
                {
                    var stringBuilder = new StringBuilder();
                    stringBuilder.Append(r[0]);
                    for (int i = 1; i < r.Length; i++)
                    {
                        stringBuilder.Append(" " + r[i]);
                    }

                    return stringBuilder.ToString();
                }).ToArray();

            await File.WriteAllLinesAsync(filename, lines);
        }
    }
}
