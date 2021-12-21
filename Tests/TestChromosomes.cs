// Created by Dr. Adam Streck, 2021, adam.streck@gmail.com

using System;
using NUnit.Framework;
using SimChA.DataTypes;

namespace Tests;

[TestFixture]
public class TestChromosomes
{
    private Chromosome _chr;
    
    [SetUp]
    public void Setup()
    {
        _chr = new Chromosome(ReferenceGenome.GetGenotype(true)[0]);
    }

    [Test]
    public void TestSplit()
    {
        Console.WriteLine(_chr.Length());
        _chr.Split(1000, true);
        Assert.AreEqual(1000, _chr.Length());
    }
    

    // public Chromosome Split(int pos, bool keepFirst)


    // public void Join(Chromosome other, bool prepend)


    // public void InvertRange(int invStart, int invEnd)

    
    // public void DuplicateRange(int start, int end)


    // public void Bridge(int pos, bool cutFront)


    // public void ScatterAndGather(List<int> locs, int count)

}