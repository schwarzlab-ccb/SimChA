// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes;

public struct Region
{
    public ChromID ChromId;
    public int Start; // Zero-index of first base
    public int End; // Zero-index one beding the last base
    public bool Forward; // True if the region is placed in the forward direction, false otherwises

    public int Length => Start - End;

    public Region(int start, int end, ChromID chromId, bool forward = true)
    {
        Start = start;
        End = end;
        ChromId = chromId;
        Forward = forward;
    }
}