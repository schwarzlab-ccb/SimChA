#!/usr/bin/env python 

import os
import argparse
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from PIL import Image
from matplotlib import cm
import random
import sys
import cv2

eps = sys.float_info.epsilon


# Take from https://stackoverflow.com/questions/3160699/python-progress-bar
def progressbar(it, prefix="", size=60, file=sys.stdout):
    count = len(it)

    def show(j):
        x = int(size * j / count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#" * x, "." * (size - x), j, count))
        file.flush()

    show(0)
    for i, item in enumerate(it):
        yield item
        show(i + 1)
    file.write("\n")
    file.flush()


def build_tree(parent_df, ids, rootId):
    tree = {cloneId: list(children) for cloneId in ids if
            len(children := parent_df.loc[parent_df["ParentId"] == cloneId, "ChildId"]) > 0}
    if -1 not in tree:
        tree[-1] = [rootId]
    return tree


def line_from_data(px_row_orig, width):
    result = []
    pop_total = sum(row_it[0] for row_it in px_row_orig)
    px_row = [[px_row_orig[i][0] * width / pop_total, px_row_orig[i][1]] for i in range(len(px_row_orig))]
    contr = 1.0
    index = 0
    new_px = np.zeros(4)
    while index < len(px_row):
        if px_row[index][0] >= eps:
            factor = contr if px_row[index][0] >= contr else px_row[index][0]
            px_row[index][0] -= factor
            contr -= factor
            new_px += px_row[index][1] * factor
            if contr < eps:
                contr = 1.0
                result.append(new_px)
                new_px = np.zeros(4)
        else:
            index += 1
    if len(result) < width:
        result.append(new_px + contr * np.ones(4) * .5)  # Add gray to the last pixel if needed
    return result


def get_pop_dict(pops_df, gen, max_pop, normalize):
    test_pop = pops_df[pops_df["Gen"] == gen]
    pop_dict = {test_pop["Id"].iat[i]: test_pop["Pop"].iat[i] for i in range(len(test_pop))}
    if not normalize:
        pop_count = pops_df[pops_df["Gen"] == gen]["Pop"].sum()
        remainder = max_pop - pop_count
        pop_dict[-1] = remainder
    return pop_dict


def rec_descend(tree, pop_dict, cloneId):
    res = []
    if cloneId in tree:
        if cloneId in pop_dict:
            pop = pop_dict[cloneId] / 2
            res += [(pop, cloneId)]
            for childId in tree[cloneId]:
                res += rec_descend(tree, pop_dict, childId)
            res += [(pop, cloneId)]
        else:
            for childId in tree[cloneId]:
                res += rec_descend(tree, pop_dict, childId)
    elif cloneId in pop_dict:
        return [(pop_dict[cloneId], cloneId)]
    return res


def get_image_pixels(tree, rootId, pop_df, gen_count, col_map, breath, max_pop, normalize=True):
    img_pixels = []
    for i in progressbar(range(gen_count), "Converting generation: "):
        pop_dict = get_pop_dict(pop_df, i, max_pop, normalize)
        if len(pop_dict) > 0:
            id_rows = rec_descend(tree, pop_dict, rootId if normalize else -1)
            pixel_row = [(id_row[0], col_map[id_row[1]]) for id_row in id_rows]
        else:
            pixel_row = [(-1, np.ones(4) * 127)]
        line_px = line_from_data(pixel_row, breath)
        img_pixels.append(line_px)
    return img_pixels


def pixels_to_img(img_pixels, res):
    px_array = np.array(img_pixels, dtype=np.uint8)
    resized = cv2.resize(px_array, [res[1], res[0]])  # Switch coordintaes as it's before the rotation
    new_image = Image.fromarray(resized).transpose(Image.ROTATE_270).transpose(Image.FLIP_LEFT_RIGHT)
    return new_image


def plot_fish(populations_df, parent_df, res, scale_up=1):
    pop_sizes = populations_df[["Gen", "Pop"]].groupby("Gen", as_index=False).sum()
    gen_count = pop_sizes["Gen"].max() + 1
    max_pop = pop_sizes["Pop"].max()
    print(f"Gens: {gen_count}, Max Pop {max_pop}")

    parents = parent_df["ParentId"].unique()
    children = parent_df["ChildId"].unique()
    root_list = np.setdiff1d(parents, children)
    if len(root_list) != 1:
        raise Exception("Failed to determine root. There must be exactly one node with no parent in the parent tree.")
    rootId = root_list[0]
    ids = np.concatenate([[rootId], children])
    tree = build_tree(parent_df, ids, rootId)

    # Create Colors
    cols = cm.rainbow(np.linspace(0, 1, len(ids))).tolist()
    random.shuffle(cols)
    col_map = {ids[i]: np.array(list(map(lambda x: (int)(x * 255), cols[i]))) for i in range(len(ids))}
    col_map[-1] = np.ones(4)  # Root is white

    img_pixels = get_image_pixels(tree, rootId, populations_df, gen_count, col_map, res[1] * scale_up, max_pop, True)
    image = pixels_to_img(img_pixels, res)
    return image


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TODO pyfish')
    parser.add_argument("populations", type=str, help="TODO pop")
    parser.add_argument("parent_tree", type=str, help="TODO tree")
    parser.add_argument("output", type=str, help="TODO output")
    args = parser.parse_args()

    res = [2160, 1080]
    random.seed(3)

    populations_df = pd.read_csv(args.populations)
    parent_df = pd.read_csv(args.parent_tree)
    img = plot_fish(populations_df, parent_df, res)
    img.save(args.output)
