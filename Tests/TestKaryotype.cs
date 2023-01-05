using System;
using NUnit.Framework;
using SimChA.DataTypes;
using SimChA.IO;
using SimChA.Simulation;

namespace Tests;

[TestFixture]
public class TestKaryotype
{
    private Karyotype _kar;
    private Random _rnd;
    
    [SetUp]
    public void Setup()
    {
        _kar = new Karyotype(false);
        _rnd = new Random(0);
    }

    [Test]
    public void TestDeletion()
    {
        _kar.ApplyAberration(_rnd, AberrationEnum.ChromDeletion, new BaseAbbP(1f));
        Assert.AreEqual(45, _kar.ChrCount);
    }
    
    [Test]
    public void TestDuplication()
    {
        _kar.ApplyAberration(_rnd, AberrationEnum.ChromDuplication, new BaseAbbP(1f));
        Assert.AreEqual(47, _kar.ChrCount);
    }
    
    [Test]
    public void TestBFB()
    {
        _kar.ApplyAberration(_rnd, AberrationEnum.BreakageFusionBridge, new FractionAbbP(1f, .1f));
        Assert.AreEqual(46, _kar.ChrCount);
    }
    
    [Test]
    public void TestClean()
    {
        for (int i = 0; i < 46; i++)
        {
            _kar.ApplyAberration(_rnd, AberrationEnum.ChromDeletion, new BaseAbbP(1f));
        }
        Assert.AreEqual("[]", _kar.ToString());
    }
}