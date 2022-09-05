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
    
    [SetUp]
    public void Setup()
    {
        _kar = new Karyotype(false, new Random(0));
    }

    [Test]
    public void TestDeletion()
    {
        _kar.ApplyAberration(AberrationEnum.ChromDeletion, new BaseAbbP(1f));
        Assert.AreEqual(45, _kar.ChromCount);
    }
    
    [Test]
    public void TestDuplication()
    {
        _kar.ApplyAberration(AberrationEnum.ChromDuplication, new BaseAbbP(1f));
        Assert.AreEqual(47, _kar.ChromCount);
    }
    
    [Test]
    public void TestBFB()
    {
        _kar.ApplyAberration(AberrationEnum.BreakageFusionBridge, new FractionAbbP(1f, .1f));
        Assert.AreEqual(46, _kar.ChromCount);
    }
    
    [Test]
    public void TestClean()
    {
        for (int i = 0; i < 46; i++)
        {
            _kar.ApplyAberration(AberrationEnum.ChromDeletion, new BaseAbbP(1f));
        }
        Assert.AreEqual("[]", _kar.ToString());
    }
}