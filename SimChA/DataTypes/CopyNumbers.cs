namespace SimChA.DataTypes;

public class CopyNumbers
{
    public List<Region> AllRegions;
    public bool IsFemale;
    public IEnumerable<ChromNum> ReferenceChromosomes;
    private List<Region> _segmentation = new List<Region>();
    private List<(int, int)> _copynumbers = new List<(int, int)>();

    public CopyNumbers(Karyotype karyotype)
    {
        AllRegions = karyotype.GetAllRegions();
        IsFemale = karyotype.IsFemale;
        ReferenceChromosomes = ReferenceGenome.GetChromosomes(IsFemale);

        // minimum consistent segmentation
        CalcSegmentation();

        // Counting
        CalcCopyNumbers();

        // TODO: Merge neighboring segments that have the same copy numbers
    }

    private void CalcSegmentation()
    {
        foreach (var curRefChrom in ReferenceChromosomes)
        {
            var curRegions = AllRegions.Where(region => region.ChromId.ChromNum == curRefChrom).ToList();
            curRegions.Add(ReferenceGenome.GetRegion(curRefChrom));

            var segmentBoundaries = curRegions.Select(r => r.Start).Concat(curRegions.Select(r => r.End)).Distinct().ToList();
            segmentBoundaries.Sort();

            for (int i = 0; i < segmentBoundaries.Count() - 1; i++)
            {
                _segmentation.Add(new Region(segmentBoundaries[i], segmentBoundaries[i + 1], new ChromID(curRefChrom, true)));
            }

        }
    }

    private void CalcCopyNumbers()
    {
        foreach (var segment in _segmentation)
        {
            var curRegions = AllRegions.Where(region => region.ChromId.ChromNum == segment.ChromId.ChromNum);
            var copynumberH1 = curRegions.Where(r => r.ChromId.Parent && segment.IsInside(r)).Count();
            var copynumberH2 = curRegions.Where(r => !r.ChromId.Parent && segment.IsInside(r)).Count();
            _copynumbers.Add((copynumberH1, copynumberH2));

            // Console.WriteLine($"{segment}");
            // foreach (var r in curRegions)
            // {
            //     Console.WriteLine($"\t{r}");
            // }
            // Console.WriteLine($"\t{curRegions.Count()} - {copynumberH1} / {copynumberH2}");
        }
    }

    public override string ToString()
    {
        var outputString = "";
        for (int i = 0; i < _segmentation.Count(); i++)
        {
            outputString = outputString + _segmentation[i] + ": " + _copynumbers[i] + "\n";
        }
        return outputString;
    }


    public string ToTSV(bool isFirst = true)
    {
        var outputString = isFirst ? "chrom\tstart\tend\tcn_a\tcn_b\n" : "";
        for (int i = 0; i < _segmentation.Count(); i++)
        {
            string[] curLine = { _segmentation[i].ChromId.ChromNum.ToString(), _segmentation[i].Start.ToString(), _segmentation[i].End.ToString(), _copynumbers[i].Item1.ToString(), _copynumbers[i].Item2.ToString() };
            outputString = outputString + string.Join('\t', curLine) + "\n";
        }
        return outputString;
    }

    public string ToTSV(string sampleId, bool isFirst = true)
    {
        var outputString = isFirst ? "sample_id\tchrom\tstart\tend\tcn_a\tcn_b\n" : "";
        for (int i = 0; i < _segmentation.Count(); i++)
        {
            string[] curLine = { sampleId, _segmentation[i].ChromId.ChromNum.ToString(), _segmentation[i].Start.ToString(), _segmentation[i].End.ToString(), _copynumbers[i].Item1.ToString(), _copynumbers[i].Item2.ToString() };
            outputString = outputString + string.Join('\t', curLine) + "\n";
        }
        return outputString;
    }

}