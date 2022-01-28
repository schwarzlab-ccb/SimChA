using System;
using NUnit.Framework;
using SimChA.DataTypes;
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
    public void TestMissegration()
    {
        var res = _kar.ApplyAbberation(AberrationEnum.Missegregation);
        Assert.AreEqual(47, _kar.ChromCount);
        Assert.AreEqual(45, res.ChromCount);
    }
    
    [Test]
    public void TestDuplication()
    {
        var res = _kar.ApplyAbberation(AberrationEnum.Duplication);
        Assert.AreEqual(47, _kar.ChromCount);
        Assert.AreEqual(46, res.ChromCount);
    }
    
    [Test]
    public void TestBFB()
    {
        var res = _kar.ApplyAbberation(AberrationEnum.BreakageFusionBridge);
        Assert.AreEqual(46, _kar.ChromCount);
        Assert.AreEqual(45, res.ChromCount);
    }
    
    [Test]
    public void TestClean()
    {
        for (int i = 0; i < 46; i++)
        {
            _kar = _kar.ApplyAbberation(AberrationEnum.Missegregation);
        }
        Assert.AreEqual("[]", _kar.ToString());
    }
}