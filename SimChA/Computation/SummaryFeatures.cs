using SimChA.DataTypes;
using SimChA.Simulation;
namespace SimChA.Computation;

public static class SummaryFeatures
{
    // Get the copy-number profiles of the clones in the sample
    public static Dictionary<Karyotype, List<CopyNumber>> GetCopyNumberProfiles(GenRef genRef, IList<Karyotype> kars)
    {
        var cnProfiles = new Dictionary<Karyotype, List<CopyNumber>>();
        foreach (var kar in kars)
        {
            var cnList = CopyNumbers.CalcCopyNumbers(genRef, kar, kar.SexXX).ToList();
            cnProfiles.Add(kar, cnList);
        }
        return cnProfiles;
    }

    public static List<long> GetSegLengths(GenRef genRef, IList<Karyotype> kars, bool includeCNNormal = false, bool includeLOH = false, bool includeSexChromosomes = false)
    {
        var segLengths = new List<long>();
        foreach (var cnProfile in GetCopyNumberProfiles(genRef, kars))
        {
            var kar = cnProfile.Key;
            var cnList = cnProfile.Value;
            // Drop the sex chromosomes
            if (!includeSexChromosomes)
            {
                cnList = cnList.Where(cn => !(cn.Segment.ChrNo == "chrX" || cn.Segment.ChrNo == "chrY")).ToList();
            }
            if (!includeCNNormal)
            {
                cnList = cnList.Where(cn => !(cn.CNH1 == 1 && cn.CNH2 == 1)).ToList();
            }
            if (!includeLOH)
            {
                cnList = cnList.Where(cn => !(cn.CNH1 + cn.CNH2 == 2 && (cn.CNH1 == 0 || cn.CNH2 == 0))).ToList();
            }
            segLengths.AddRange(cnList.Select(cn => cn.Segment.Length));
        }
        return segLengths;
    }

    public static Dictionary<string, List<long>> GetChangepoints(GenRef genRef, IList<Karyotype> kars)
    {
        var changepointDict = genRef.AllChrs.ToDictionary(chr => chr, chr => new List<long>());
        return changepointDict;
    }
}