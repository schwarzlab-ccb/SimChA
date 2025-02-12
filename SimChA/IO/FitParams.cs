// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

namespace SimChA.IO;

public record FitParams
(
    double Stress = 1,
    double TsgOg = 1,
    double Essentiality = 1,
    bool GeneNormalization = false
);
