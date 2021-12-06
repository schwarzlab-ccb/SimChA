// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

namespace SimChA.DataTypes
{
    public struct Region
    {
        public int Start; // Zero-index of first base
        public int End; // Zero-index one beding the last base
        public ChromID ChromId;

        public int Length => Start - End;

        public Region(int start, int end, ChromID chromId)
        {
            Start = start;
            End = end;
            ChromId = chromId;
        }
    }
}