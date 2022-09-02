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
    public void TestDeletion()
    {
        _kar.ApplyAbberation(AberrationEnum.ChromDeletion);
        Assert.AreEqual(45, _kar.ChromCount);
    }
    
    [Test]
    public void TestDuplication()
    {
        _kar.ApplyAbberation(AberrationEnum.ChromDuplication);
        Assert.AreEqual(47, _kar.ChromCount);
    }
    
    [Test]
    public void TestBFB()
    {
        _kar.ApplyAbberation(AberrationEnum.BreakageFusionBridge);
        Assert.AreEqual(46, _kar.ChromCount);
    }
    
    [Test]
    public void TestClean()
    {
        for (int i = 0; i < 46; i++)
        {
            _kar.ApplyAbberation(AberrationEnum.ChromDeletion);
        }
        Assert.AreEqual("[]", _kar.ToString());
    }
}