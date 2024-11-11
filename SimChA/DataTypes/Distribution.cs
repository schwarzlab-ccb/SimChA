// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System.Text.Json.Serialization;

namespace SimChA.DataTypes;

[JsonConverter(typeof(JsonStringEnumConverter))]
public enum Distribution { Uniform, Exponential, Normal, Unit, Poisson, Geometric }
