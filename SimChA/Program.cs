using MathNet.Numerics.Distributions;
using SimChA.DataTypes;

Console.WriteLine("SimChA 0.0.1");
var cloneList = new CloneList();
Console.WriteLine(cloneList);

int binomial = Binomial.Sample(0.5, 100);
Console.WriteLine(binomial);