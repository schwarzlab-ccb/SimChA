// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

namespace SimChA.DataTypes;

public class GenRef
{
    public string Name { get; }
    public Dictionary<ChrNo, int> ChrLengths { get; }
    public Dictionary<ChrNo, SexEnum> ChrSex { get; }
    public int AutosomeCount { get; }
    public int ChrCount { get; }
    
    private Region[] XYGenome { get; }
    private Region[] XXGenome { get; }
    private long XYGenomeLen { get; }
    private long XXGenomeLen { get; }
    
    public Region[] GetGenotype(bool sexXX)
        => sexXX ? XXGenome : XYGenome;

    public long GetGenomeLen(bool sexXX)
        => sexXX ? XXGenomeLen : XYGenomeLen;

    public long GetChromLen(ChrNo chrNo)
        => ChrLengths[chrNo];
    
    public IEnumerable<ChrNo> ChrIDsForSex(bool sexXX)
        => Enum.GetValues<ChrNo>().Take(sexXX ? AutosomeCount + 1 : AutosomeCount + 2);
    
    public GenRef(string name, Dictionary<ChrNo, int> chrLengths, Dictionary<ChrNo, SexEnum> chrSex)
    {
        Name = name;
        ChrLengths = chrLengths;
        ChrSex = chrSex;
        AutosomeCount = chrSex.Count(x => x.Value == SexEnum.Both);
        ChrCount = AutosomeCount * 2 + (chrSex.Count - AutosomeCount);
        var haplotypeOneF = CreateHaplotype(true, true);
        var haplotypeTwoF = CreateHaplotype(false, true);
        var haplotypeOneM = CreateHaplotype(true, false);
        var haplotypeTwoM = CreateHaplotype(false, false);
        XYGenome = haplotypeOneM.Concat(haplotypeTwoM).ToArray();
        XXGenome = haplotypeOneF.Concat(haplotypeTwoF).ToArray();
        XYGenomeLen = XYGenome.Sum(r => r.Length);
        XXGenomeLen = XXGenome.Sum(r => r.Length);
    }
    
    private Region GetRegion(ChrNo chrNo, bool isFirstHaplotype = true) 
        => new(0, ChrLengths[chrNo], chrNo, isFirstHaplotype);
    
    private IEnumerable<Region> CreateHaplotype(bool isFirstHaplotype, bool isFemale)
    {
        var nonGender = ChrSex.Select(x => x.Key).Where(x => ChrSex[x] == SexEnum.Both);
        var sexChr = isFirstHaplotype | isFemale ? ChrNo.chrX : ChrNo.chrY;
        var all = nonGender.Concat(new[] { sexChr });
        return all.Select(num => GetRegion(num, isFirstHaplotype));
    }
}