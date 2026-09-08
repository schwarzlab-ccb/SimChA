// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using SimChA.Computation;

namespace SimChA.IO;

public record FitParams
(
    double Stress = 0,
    double TsgOg = 0,
    double Essentiality = 0,
    // Shapes the essentiality penalty across copy number: what losing one of two copies costs
    // relative to losing both. See Fitness.DosageLoss -- 1 is linear, 2 charges a quarter for a
    // single copy, and a large value reproduces the CN==0-only step used before this was a
    // parameter. Configs written before it existed deserialise to this default.
    double HaploExponent = 2.0,
    string GeneSet = "Empty"
);
