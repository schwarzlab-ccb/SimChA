// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public struct ChromID
{
    public ChromNum ChromNum;
    public bool Parent;

    public ChromID(ChromNum num, bool parent)
    {
        ChromNum = num;
        Parent = parent;
    }

    public override string ToString()
    {
        string parentIdentifier = Parent ? "M" : "P";
        return $"({ChromNum.ToString()},{parentIdentifier})";
    }
}