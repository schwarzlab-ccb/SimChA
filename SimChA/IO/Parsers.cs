using System.Collections;
using System.Globalization;
using System.Text.Json;
using SimChA.DataTypes;
using SimChA.Simulation;
using System.Text;
using System.Text.RegularExpressions;

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

    // Expected format is that there is a header and the columns contain:
    // SampleID, Chr, Start, End, CN hap1, CN hap2
    // NOTE: This became quite unwieldy due to the missing regions calculation,
    // however it works so don't refactor unless needed
    public static Dictionary<string, Karyotype> ParseCNAProfile(GenRef genRef, TextReader cnaFile)
    {
        Dictionary<string, Karyotype> result = new();
        var missingRanges = new List<GenRange>();
        var regionsA = new List<Region>();
        var regionsB = new List<Region>();
        string lastSample = "";
        var lastChr = ChrNo.chr1;
        long lastPos = 0L;
        cnaFile.ReadLine(); // Skip header
        var present = Enum.GetValues<ChrNo>().ToDictionary(c => c, _ => false);
        while (cnaFile.ReadLine() is { } line)
        {
            string[] lineSplit = line.Split('\t');
            string sample = lineSplit[0];
            // Set the new sample
            if (sample != lastSample)
            {
                // First is empty
                if (regionsA.Any() || regionsB.Any())
                {
                    // Till the end of a chromosome
                    if (lastPos != genRef.ChrLengths[lastChr])
                    {
                        missingRanges.Add(new GenRange(lastPos, genRef.ChrLengths[lastChr], lastChr));
                    }
                    missingRanges.AddRange(
                        present
                            .Where(pair => !pair.Value)
                            .Select(c => new GenRange(0, genRef.ChrLengths[c.Key], c.Key)));
                    // Consider missing to be haplotypes by default
                    foreach (var range in missingRanges)
                    {
                        regionsA.Add(new Region(range.Start, range.End, range.ChrNo, true));
                        regionsB.Add(new Region(range.Start, range.End, range.ChrNo, false));
                    }
                    bool sexXX = !present[ChrNo.chrY];
                    result[lastSample] = new Karyotype(new List<Contig> { new(regionsA),  new(regionsB) }, missingRanges, sexXX);
                }
                // Reset
                regionsA.Clear();
                regionsB.Clear();
                missingRanges.Clear();
                lastSample = sample;
                lastChr = ChrNo.chr1;
                lastPos = 0L;
                foreach (var pair in present)
                {
                    present[pair.Key] = false;
                }
            }
            try
            {
                // Parse the line
                var chrNo = (ChrNo)Enum.Parse(typeof(ChrNo), lineSplit[1]);
                present[chrNo] = true;
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
                    if (lastPos != genRef.ChrLengths[lastChr])
                    {
                        missingRanges.Add(new GenRange(lastPos, genRef.ChrLengths[lastChr], lastChr));
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
                    regionsA.Add(new Region(start, end, chrNo, true));
                }
                for (int i = 0; i < cnB; i++)
                {
                    regionsB.Add(new Region(start, end, chrNo, false));
                }
            }
            catch (Exception e)
            {
                throw new Exception($"Could not parse the CNA profile:\n{line}\n{e.Message}");
            }
        }

        // Till the end of a chromosome
        if (lastPos != genRef.ChrLengths[lastChr])
        {
            missingRanges.Add(new GenRange(lastPos, genRef.ChrLengths[lastChr], lastChr));
        }
        missingRanges.AddRange(
            present
                .Where(pair => !pair.Value)
                .Select(c => new GenRange(0, genRef.ChrLengths[lastChr], c.Key)));
        // Consider missing to be haplotypes by default
        foreach (var range in missingRanges)
        {
            regionsA.Add(new Region(range.Start, range.End, range.ChrNo, true));
            regionsB.Add(new Region(range.Start, range.End, range.ChrNo, false));
        }

        // Add the last sample
        bool sexXX1 = !present[ChrNo.chrY];
        result[lastSample] = new(new List<Contig> { new(regionsA),  new(regionsB) }, missingRanges, sexXX1);
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
            var chrNum = (ChrNo)Enum.Parse(typeof(ChrNo), genString[0]);
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
        var columns = new Dictionary<string, int> { { idKey, -1 }, { parentIDKey, -1 }, { distanceKey, -1 } };
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

    public static (Dictionary<ChrNo, int> chrLengths, Dictionary<ChrNo, SexEnum> chrSex) ParseChromosomes(string text)
    {
        IList<string> lines = text.Split("\n");
        Dictionary<ChrNo, int> chrLengths = new();
        Dictionary<ChrNo, SexEnum> chrSex = new();
        for (var index = 0; index < lines.Count; index++)
        {
            string line = lines[index];
            var lineSplit = line.Split("\t").Select(s => s.Trim()).ToList();
            var chrNo = (ChrNo) Enum.Parse(typeof(ChrNo), lineSplit[0]);
            int length = int.Parse(lineSplit[1]);
            chrLengths.Add(chrNo, length);
            var sexEnum = GetSexEnum(lines, lineSplit, index);
            chrSex.Add(chrNo, sexEnum);
        }

        return (chrLengths, chrSex);
    }

    private static SexEnum GetSexEnum(IList<string> lines, IReadOnlyList<string> lineSplit, int index)
    {
        if (lineSplit.Count > 2)
        {
            string sexString = lineSplit[2];
            return Enum.Parse<SexEnum>(sexString);
        }
        if (index == lines.Count - 1)
        {
            return SexEnum.Male;
        }
        if (index == lines.Count - 2)
        {
            return SexEnum.Female;
        }
        return SexEnum.Both;
    }
    
    public static IEnumerable<GenContents> ParseFasta(StreamReader fastaStream)
    {
        GenContents? genContents = null;
        while (fastaStream.ReadLine() is { } line)
        {
            if (line.StartsWith(";"))
            {
                continue;
            }
            if (line.StartsWith(">"))
            {
                if (genContents != null)
                {
                    yield return genContents;
                }
                string pattern = @"^>chr([1-9]|1[0-9]|2[0-2]|X|Y)$";
                var match = Regex.Match(line, pattern);
                if (match.Value == "" || !Enum.TryParse(match.Value[1..], out ChrNo chrNo))
                {
                    genContents = null;
                    continue;
                }
                Console.WriteLine(chrNo);
                genContents = new GenContents{ChrNo = chrNo, Sequence = new StringBuilder("", (int)HGRef.GetChromLen(chrNo))};
            }
            else if (genContents != null)
            {
                genContents.Sequence.Append(line);
            }
        }
        if (genContents != null)
        {
            yield return genContents;
        }
    }
}