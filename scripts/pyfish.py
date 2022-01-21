import os
import argparse
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image


def create_plot(pop_df, adj_df):
    x_len = pop_df["Generation"].max() + 1
    y_len = 1000

    pixels = [
       [(54, 54, 54), (232, 23, 93), (71, 71, 71), (168, 167, 167)],
       [(204, 82, 122), (54, 54, 54), (168, 167, 167), (232, 23, 93)],
       [(71, 71, 71), (168, 167, 167), (54, 54, 54), (204, 82, 122)],
       [(168, 167, 167), (204, 82, 122), (232, 23, 93), (54, 54, 54)]
    ]

    # Convert the pixels into an array using numpy
    array = np.array(pixels, dtype=np.uint8)

    # Use PIL to create an image from the new array of pixels
    new_image = Image.fromarray(array)
    return new_image


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TODO pyfish')
    parser.add_argument("populations", type=str, help="TODO pop")
    parser.add_argument("adjacency", type=str, help="TODO adj")
    args = parser.parse_args()

    populations_df = pd.read_csv(args.populations)
    adjacency_df = pd.read_csv(args.adjacency)

    img = create_plot(populations_df, adjacency_df)
    plt.imshow(img)
    # plt.show()
