namespace SimChA.DataTypes;

public struct CopyNumber
{
    public Region Segment;
    public int CNH1;
    public int CNH2;

    public string ToTSV()
        => string.Join('\t', Segment.ChromId.ChromNum.ToString(), Segment.Start.ToString(), Segment.End.ToString(), CNH1.ToString(), CNH2.ToString());
}