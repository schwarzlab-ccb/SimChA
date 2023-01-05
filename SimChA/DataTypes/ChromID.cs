// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public record ChromID(ChromNum ChromNum, bool Parent)
{
    public override string ToString()
    {
        string parentIdentifier = Parent ? "H1" : "H2";
        return $"({ChromNum.ToString()},{parentIdentifier})";
    }
}