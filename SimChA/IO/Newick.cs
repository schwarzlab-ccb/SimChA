using SimChA.DataTypes;
using SimChA.Simulation;
using System.Text.RegularExpressions;

namespace SimChA.IO;

public static class Newick
{
    private static Clone CreateNodes(string newickNode, int parentId, bool isFemale)
    {
        string[] cloneString = newickNode.Split(':');
        var clone = new Clone(int.Parse(cloneString[0].Split('-')[0]), parentId, "1", int.Parse(cloneString[1]), new Karyotype(isFemale));
        return clone;
    }


    public static List<Clone> ParseNewick(string newickString, bool isFemale)
    {
        List<Clone> clones = new();
        string regexPattern = @"(?<closeChildren>[(])|" +
                              @"(?<openChildren>[)])|" +
                              @"(?<branchLength>:+[0-9a-zA-Z-_.]+)|" +
                              @"(?<nodeName>[0-9a-zA-Z-_.]+)|" +
                              @"(?<nextNode>[,])|" +
                              @"(?<root>[;])";
        //Reverse order of newick file to start with root
        var regexOptions = RegexOptions.RightToLeft | RegexOptions.IgnorePatternWhitespace;
        if(newickString == "")
        {
            return clones;
        }
        var matches = Regex.Matches(newickString, regexPattern, regexOptions);
        if(!matches.Any() || matches[0].Value != ";") 
        {
            throw new Exception("Newick file is not in the right format");
        }
        var branchLength = CheckBranchLength(matches);
        //Iterate throu Regex-Matches
        var parentIds = new List<int> {-1};
        foreach(Match match in matches)
        {
            
            switch(match.Value)
            {
                case ";":
                    //create root node
                    clones.Add(CreateClone(clones.Count(), parentIds.Last(), match.NextMatch(), 
                        match.NextMatch().NextMatch(), isFemale, branchLength));
                    parentIds.Add(clones[clones.Count()-1].CloneId);
                    break;
                case ")":
                    //create new children
                    parentIds.Add(clones[clones.Count()-1].CloneId);
                    clones.Add(CreateClone(clones.Count(), parentIds.Last(), match.NextMatch(), 
                        match.NextMatch().NextMatch(), isFemale, branchLength));
                    clones[parentIds.Last()].ChildrenIDs.Add(clones.Count()-1);
                    break;
                case "(":
                    //remove parent from parentList and add childrenID to parent
                    parentIds.RemoveAt(parentIds.Count()-1);
                    break;
                case ",":
                    //create new child
                    clones.Add(CreateClone(clones.Count(), parentIds.Last(), match.NextMatch(), 
                        match.NextMatch().NextMatch(), isFemale, branchLength));
                    clones[parentIds.Last()].ChildrenIDs.Add(clones.Count()-1);
                    break;
            }
        }
        if(!clones.Any())
        {
            throw new Exception("No clones found in newick file. Might not be the right format.");
        }
        return clones;
    }
    
    //create clone from newick Match
    private static Clone CreateClone(int id, int parentId, Match branchLengthMatch, Match nameMatch, bool isFemale, bool branchLength)
    {
        string nameClone = nameMatch.Groups["nodeName"].Value != "" ? nameMatch.Value : 
            branchLengthMatch.Groups["nodeName"].Value != "" ? branchLengthMatch.Value : "C" + id.ToString();
        int mutCount = branchLengthMatch.Groups["branchLength"].Value != "" ? 
            (int)Math.Ceiling(float.Parse(branchLengthMatch.Value.Remove(0,1))) : branchLength ? 0 : 1;
        var clone = new Clone(id, parentId, nameClone, mutCount, new Karyotype(isFemale));
        return clone;
    }
    
    //check for branch-length in newick file
    private static bool CheckBranchLength(MatchCollection matches)
    {
        bool branchLength = false;
        foreach (Match match in matches)
        {
            if(match.Groups["branchLength"].Value != "")
            {
                branchLength = true;
                break;
            }
        }
        Console.Write(branchLength ? "" : "No branch-lengths were found, using 1 as branch-length.");
        return branchLength;
    }
}