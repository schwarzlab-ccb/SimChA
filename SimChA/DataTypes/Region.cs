// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public struct Region
{
    public ChromID ChromId;
    public int Start; // Zero-based-index of first position
    public int End; // Zero-based-index of one beyond the last position
    public bool Forward; // True if the region is placed in the forward direction, false otherwise

    public int Length => End - Start;

    public Region(int start, int end, ChromID chromId, bool forward = true)
    {
        Start = start;
        End = end;
        ChromId = chromId;
        Forward = forward;
    }

    public bool IsInside(Region other) => Start >= other.Start && End <= other.End;

    private string DirString => Forward ? "+" : "-";

    public override string ToString() => $"{ChromId}{DirString}[{Start}:{End})";
}