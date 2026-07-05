using SimChA.Data;

namespace SimChA.Computation;

public static class CopyNumbers
{
    public static Dictionary<string, List<int>> GetJointSegmentation(List<string> chromNames, List<Sample> samples)
    {
        var breaks = samples.Select(s => s.Karyotype.CalcBreaks()).ToList();
        var segmentation = chromNames.ToDictionary(
            chrom => chrom, 
            chrom => breaks.SelectMany(br => br.GetValueOrDefault(chrom, [])).ToHashSet().OrderBy(b => b).ToList());
        return segmentation;
    }

    public static List<CopyNumber> CalcCNs(IDictionary<string, List<int>> allBreaks, List<Contig> contigs)
    {
        var result = new List<CopyNumber>();
        foreach ((string chrom, var breaks) in allBreaks)
        {
            for (int i = 0; i < breaks.Count - 1; i++)
            {
                int start = breaks[i];
                int end = breaks[i + 1];
                var seg = new GenRange(start, end, chrom);
                var cns = contigs.Select(c => c.GetCNs(seg));
                (int cnA, int cnB, int nSNVs) = cns.Aggregate((CNA: 0, CNB: 0, SNV: 0), (acc, vals)
                    => (acc.CNA + vals.CNA, acc.CNB + vals.CNB, acc.SNV + vals.SNV));
                var cn = new CopyNumber(start, end, chrom, cnA, cnB, nSNVs);
                result.Add(cn);
            }
        }
        return result;
    }

    public static List<CopyNumber> CalcCNs(Karyotype karyotype, IDictionary<string, List<int>>? breaks = null)
        => karyotype.CalcCNs(breaks ?? karyotype.CalcBreaks());

    // Reports the per-base copy-number change between two karyotypes as signed deltas over
    // reference coordinates: each returned segment carries the H1/H2 change (after - before),
    // adjacent segments with the same change are merged, and unchanged segments are dropped.
    // Because it compares copy number per reference position rather than region identities, an
    // event that splits a region (e.g. an internal deletion) reports only the bases actually
    // lost or gained instead of a spurious loss plus two gains.
    public static List<CopyNumber> DiffKaryotypes(Karyotype before, Karyotype after)
    {
        var beforeBreaks = before.CalcBreaks();
        var afterBreaks = after.CalcBreaks();
        // Union both karyotypes' breakpoints so no segment straddles a change on either side.
        var jointBreaks = beforeBreaks.Keys.Union(afterBreaks.Keys).ToDictionary(
            chrom => chrom,
            chrom => beforeBreaks.GetValueOrDefault(chrom, [])
                .Concat(afterBreaks.GetValueOrDefault(chrom, []))
                .ToHashSet().OrderBy(b => b).ToList());

        // Both use the same jointBreaks instance, so the segment lists align index-for-index.
        var beforeCNs = CalcCNs(before, jointBreaks);
        var afterCNs = CalcCNs(after, jointBreaks);

        var deltas = new List<CopyNumber>();
        CopyNumber? run = null;
        for (int i = 0; i < beforeCNs.Count; i++)
        {
            var seg = beforeCNs[i];
            int dH1 = afterCNs[i].CNH1 - seg.CNH1;
            int dH2 = afterCNs[i].CNH2 - seg.CNH2;
            if (run is not null && run.Chrom == seg.Chrom && run.End == seg.Start
                && run.CNH1 == dH1 && run.CNH2 == dH2)
            {
                run = new CopyNumber(run.Start, seg.End, run.Chrom, dH1, dH2, 0);
            }
            else
            {
                if (run is { } r && (r.CNH1 != 0 || r.CNH2 != 0)) deltas.Add(r);
                run = new CopyNumber(seg.Start, seg.End, seg.Chrom, dH1, dH2, 0);
            }
        }
        if (run is { } last && (last.CNH1 != 0 || last.CNH2 != 0)) deltas.Add(last);
        return deltas;
    }

    public static double CalcPloidy(RefGen refGen, List<CopyNumber> copyNumbers, SexType sex)
        => 2 * copyNumbers.Select(c => c.Length * (c.CNH1 + c.CNH2)).Sum() / (float) refGen.SexGenomeLen[(int)sex];
}
