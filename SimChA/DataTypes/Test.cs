class A
{
    public A(int param)
    {
        // Constructor logic for A
    }
}

class B : A
{
    public B(int param) : base(param)
    {
        // Constructor logic for B
    }
}

class myClass<T> where T : B, new()
{
    public T GetElem(int param)
    {
        return new B(param);
    }
}