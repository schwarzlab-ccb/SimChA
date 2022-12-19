namespace SimChA.DataTypes;

public struct CopyNumber
{
    public Region Segment;
    public int CNH1;
    public int CNH2;

    public string ToTSV()
        => string.Join('\t', Segment.ChromId.ChromNum, Segment.Start, Segment.End, CNH1, CNH2);
}