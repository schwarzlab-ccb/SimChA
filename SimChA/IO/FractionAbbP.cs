// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.IO;

[Serializable]
public record FractionAbbP(double Likelihood, double MeanLength) : BaseAbbP(Likelihood);