using System.Globalization;
using System.Text.Json;
using SimChA.Computation;
using SimChA.DataTypes;
using SimChA.Simulation;

namespace SimChA.IO;

public static class Parsers
{
    public static SimParams ParseSimParams(string serializedJSON)
    {
        SimParams? res;
        var options = new JsonSerializerOptions { IncludeFields = true };
        try
        {
            res = JsonSerializer.Deserialize<SimParams>(serializedJSON, options);
            if (res is null)
            {
                throw new Exception($"Could not parse the simulation parameters:\n{serializedJSON}");
            }
        }
        catch (JsonException)
        {
            throw new Exception($"Could not parse the simulation parameters:\n{serializedJSON}");
        }
        if (res.Seed < 0)
        {
            return res with { Seed = new Random().Next() };
        }
        return res;
    }

    private static Karyotype MakeKaryotype(List<Region> regionsA, List<Region> regionsB, IEnumerable<GenRange> missingRanges, bool sexXX)
    {
        regionsA = RegionOps.StitchRegions(regionsA);
        var contigA = new Contig(regionsA);
        regionsB = RegionOps.StitchRegions(regionsB);
        var contigB = new Contig(regionsB);
        // Add full missing chromosomes
        var chrPresent = HGRef.ChrIDsForSex(sexXX).ToDictionary(c => c, _ => false);
        regionsA.ForEach(r => chrPresent[r.ChrNo] = true);
        regionsB.ForEach(r => chrPresent[r.ChrNo] = true);
        var missingChrs = chrPresent.Where(pair => !pair.Value).Select(pair => new GenRange(0, HGRef.GetChromLen(pair.Key), pair.Key));
        var totalMissing = missingChrs.Concat(missingRanges).ToList();
        return new Karyotype(new List<Contig> {contigA, contigB}, totalMissing, sexXX);
    }
    
    public static Dictionary<string, Karyotype> ParseCNAProfile(TextReader cnaFile)
    {
        Dictionary<string, Karyotype> result = new();
        var missingRanges = new List<GenRange>();
        var regionsA = new List<Region>();
        var regionsB = new List<Region>();
        bool sexXX = true;
        string lastSample = "";
        var lastChr = ChrNo.chr1;
        long lastPos = 0L;
        cnaFile.ReadLine(); // Skip header
        while (cnaFile.ReadLine() is { } line)
        {
            string[] lineSplit = line.Split('\t');
            string sample = lineSplit[0];
            // Set the new sample
            if (sample != lastSample)
            {
                if (regionsA.Any() || regionsB.Any())
                {
                    result[lastSample] = MakeKaryotype(regionsA, regionsB, missingRanges, sexXX);
                }
                // Reset
                regionsA.Clear();
                regionsB.Clear();
                missingRanges.Clear();
                lastSample = sample;
                lastChr = ChrNo.chr1;
                lastPos = 0L;
                sexXX = true;
            }

            // Parse the line
            var chrNo = (ChrNo) Enum.Parse(typeof(ChrNo), lineSplit[1]);
            if (chrNo == ChrNo.chrY)
            {
                sexXX = false;
            }
            int start = int.Parse(lineSplit[2]) - 1;
            int end = int.Parse(lineSplit[3]);
            int cnA = int.Parse(lineSplit[4]);
            int cnB = int.Parse(lineSplit[5]);
            
            // Check for missing ranges
            if (chrNo == lastChr)
            {
                // Range skipped
                if (lastPos != start)
                {
                    missingRanges.Add(new GenRange(lastPos, start, lastChr));
                }
            }
            else
            {
                // Till the end of a chromosome
                if (lastPos != HGRef.GetChromLen(lastChr))
                {
                    missingRanges.Add(new GenRange(lastPos, HGRef.GetChromLen(lastChr), lastChr));
                }
                // Start of a chromosome
                if (start != 0)
                {
                    missingRanges.Add(new GenRange(0, start, chrNo));
                }
            }
            lastChr = chrNo;
            lastPos = end;
            
            // Add the new regions
            for (int i = 0; i < cnA; i++)
            {
                regionsA.Add(new Region(start, end, new ChrID(chrNo, true)));
            }
            for (int i = 0; i < cnB; i++)
            {
                regionsB.Add(new Region(start, end, new ChrID(chrNo, false)));
            }
        }
        
        // Add the last sample
        result[lastSample] = MakeKaryotype(regionsA, regionsB, missingRanges, sexXX);
        return result;
    }

    public static Dictionary<ChrNo, List<Gene>> ParseGeneList(TextReader geneFile)
    {
        // Pre-initialization
        var noEnum = Enum.GetValues(typeof(ChrNo)).Cast<ChrNo>().ToList();
        var geneList = noEnum.ToDictionary(c => c, _ => new List<Gene>());
        while (geneFile.ReadLine() is { } line)
        {
            if (line == "") continue;
            string[] genString = line.Split('\t');
            string name = genString[3];
            double fitness = double.Parse(genString[4], CultureInfo.InvariantCulture.NumberFormat);
            var chrNum = (ChrNo) Enum.Parse(typeof(ChrNo), genString[0]);
            // Convert to zero-based [start, end) index 
            var region = new GenRange(int.Parse(genString[1]) - 1, int.Parse(genString[2]), chrNum);
            var gene = new Gene(name, region, fitness);
            geneList[chrNum].Add(gene);
        }
        foreach (var pair in geneList)
        {
            pair.Value.Sort((g1, g2) => g1.Range.Start.CompareTo(g2.Range.Start));
        }
        return geneList;
    }
    
    public static List<CloneIn> ParseClones(TextReader cloneStream, bool parseFitness)
    {
        const string idKey = "ID";
        const string parentIDKey = "ParentID";
        const string distanceKey = "Distance";
        const string fitnessKey = "Fitness";
        
        string? firstLine = cloneStream.ReadLine();
        if (firstLine == null) throw new Exception("CloneIn file is empty.");
        var header = firstLine.Split(",").Select(s => s.Trim()).ToList();
        var columns = new Dictionary<string, int> {{idKey, -1}, {parentIDKey, -1},  {distanceKey, -1}};
        if (parseFitness)
        {
            columns.Add(fitnessKey, -1);
        }
        
        foreach (var column in columns)
        {
            int idx = header.IndexOf(column.Key);
            if (idx == -1) throw new Exception($"CloneIn file does not contain {column.Key} column.");
            columns[column.Key] = idx;
        }

        var cloneFitness = new List<CloneIn>();
        while (cloneStream.ReadLine() is { } line)
        {
            var lineSplit = line.Split(",").Select(s => s.Trim()).ToList();
            int id = int.Parse(lineSplit[columns[idKey]]);
            int parentId = int.Parse(lineSplit[columns[parentIDKey]]);
            int distance = int.Parse(lineSplit[columns[distanceKey]]);
            double fitness = parseFitness
                ? double.Parse(lineSplit[columns[fitnessKey]], CultureInfo.InvariantCulture.NumberFormat)
                : 0.0;
            var clone = new CloneIn(id, parentId, distance, fitness);
            cloneFitness.Add(clone);
        }
        return cloneFitness;
    }
}