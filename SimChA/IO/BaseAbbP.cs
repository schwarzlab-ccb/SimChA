// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using SimChA.DataTypes;

namespace SimChA.IO;

[Serializable]
public class BaseAbbP
{
    public AberrationEnum Aberration;
    public float Likelihood;
}