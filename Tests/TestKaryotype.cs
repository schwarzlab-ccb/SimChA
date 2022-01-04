using System;
using NUnit.Framework;
using SimChA.DataTypes;

namespace Tests;

[TestFixture]
public class TestKaryotype
{
    private Karyotype _kar;
    
    [SetUp]
    public void Setup()
    {
        _kar = new Karyotype(false);
    }

    [Test]
    public void TestMissegration()
    {
        var res = _kar.ApplyAbberation(AbberationEnum.Missegregation);
        Assert.AreEqual(47, _kar.ChromCount);
        Assert.AreEqual(45, res.ChromCount);
    }

    [Test]
    public void TestClean()
    {
        for (int i = 0; i < 46; i++)
        {
            _kar = _kar.ApplyAbberation(AbberationEnum.Missegregation);
            _kar.Clean();
        }
        Console.WriteLine(_kar);
        Assert.AreEqual("[]", _kar.ToString());
    }
    
    [Test]
    public void TestDuplication()
    {
        var res = _kar.ApplyAbberation(AbberationEnum.Duplication);
        Assert.AreEqual(47, _kar.ChromCount);
        Assert.AreEqual(46, res.ChromCount);
    }
}