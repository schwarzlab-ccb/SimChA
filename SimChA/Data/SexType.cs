// Created by Dr. Adam Streck, 2023, adam.streck@gmail.com

using System.Text.Json.Serialization;

namespace SimChA.Data;

[JsonConverter(typeof(JsonStringEnumConverter))]
public enum SexType { Female, Male, None }