// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System.Text.Json.Serialization;

namespace SimChA.DataTypes;

[JsonConverter(typeof(JsonStringEnumConverter))]
public enum SexEnum
{
    Female,
    Male,
    Both
}