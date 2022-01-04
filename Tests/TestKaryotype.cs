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
        var _kar2 = _kar.ApplyAbberation(AbberationEnum.Missegregation);
        Assert.AreEqual(47, _kar.ChromCount);
        Assert.AreEqual(45, _kar2.ChromCount);
    }

    [Test]
    public void TestClean()
    {
        for (int i = 0; i < 10000; i++)
        {
            _kar.ApplyAbberation(AbberationEnum.InternalDeletion);
        }
        _kar.Clean();
        Console.Write(_kar);
    }
}