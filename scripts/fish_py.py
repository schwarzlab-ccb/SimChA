#!/usr/bin/env python

import argparse
import random
import sys
import math
import logging

import numpy as np
import pandas as pd
from PIL import Image
from matplotlib import cm

log = logging.getLogger()
log.setLevel(logging.INFO)
eps = sys.float_info.epsilon
random.seed(3)

# Taken from https://stackoverflow.com/questions/3160699/python-progress-bar
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


def build_tree(par_df, ids, rootId):
    tree = {cloneId: list(children) for cloneId in ids if
            len(children := par_df.loc[par_df["ParentId"] == cloneId, "ChildId"]) > 0}
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


# Interpolates between generations
def get_count_at_gen(pops_df, first_step, last_step, step, clone_id):
    cp_df = pops_df[pops_df["Id"] == clone_id]
    match = cp_df[cp_df["Step"] == step]
    if not match.empty:
        return float(cp_df["Pop"][match.index])

    l_neigh_i = cp_df[cp_df["Step"] < step]["Step"]
    u_neigh_i = cp_df[cp_df["Step"] > step]["Step"]

    if l_neigh_i.size > 0 and u_neigh_i.size > 0:
        ind_min = l_neigh_i.idxmax()
        ind_max = u_neigh_i.idxmin()
        prev_gen, prev_pop = cp_df.loc[ind_min, ["Step", "Pop"]].tolist()
        next_gen, next_pop = cp_df.loc[ind_max, ["Step", "Pop"]].tolist()
    elif l_neigh_i.size > 0:
        ind_min = l_neigh_i.idxmax()
        prev_gen, prev_pop = cp_df.loc[ind_min, ["Step", "Pop"]].tolist()
        next_gen, next_pop = last_step, 0
    elif u_neigh_i.size > 0:
        ind_max = u_neigh_i.idxmin()
        prev_gen, prev_pop = first_step, 0
        next_gen, next_pop = cp_df.loc[ind_max, ["Step", "Pop"]].tolist()
    else:
        return 0

    pop = ((step - prev_gen) * next_pop + (next_gen - step) * prev_pop) / (next_gen - prev_gen)
    log.debug(f"Id: {clone_id}, Step: {step}, Pop: {pop}")
    return pop


def get_pop_dict(pops_df, ids, fg, lg, gen, max_pop, normalize):
    pop_dict = { cloneId: count for cloneId in ids if (count := get_count_at_gen(pops_df, fg, lg, gen, cloneId)) > 0 }
    if not normalize:
        pop_count = sum(pop_dict.values())
        remainder = max_pop - pop_count
        pop_dict[-1] = remainder
    return pop_dict


def rec_descend(tree, pop_dict, clone_id):
    res = []
    if clone_id in tree:
        if clone_id in pop_dict:
            pop = pop_dict[clone_id] / 2
            res += [(pop, clone_id)]
            for childId in tree[clone_id]:
                res += rec_descend(tree, pop_dict, childId)
            res += [(pop, clone_id)]
        else:
            for childId in tree[clone_id]:
                res += rec_descend(tree, pop_dict, childId)
    elif clone_id in pop_dict:
        return [(pop_dict[clone_id], clone_id)]
    return res


def pixels_to_img(img_pixels, res):
    px_array = np.array(img_pixels, dtype=np.uint8)
    new_image = Image.fromarray(px_array).transpose(Image.ROTATE_270).transpose(Image.FLIP_LEFT_RIGHT)
    resized = new_image.resize(res)
    return resized


def get_image_pixels(tree, ids, root_id, pop_df, col_map, max_pop, breath, first_step, last_step, normalize=True):
    img_pixels = []
    for step in progressbar(range(first_step, last_step), "Step: "):
        pop_dict = get_pop_dict(pop_df, ids, first_step, last_step, step, max_pop, normalize)
        if len(pop_dict) > 0:
            id_rows = rec_descend(tree, pop_dict, root_id if normalize else -1)
            pixel_row = [(id_row[0], col_map[id_row[1]]) for id_row in id_rows]
        else:
            pixel_row = [(-1, np.ones(4) * 127)]
        line_px = line_from_data(pixel_row, breath)
        img_pixels.append(line_px)
    return img_pixels


def plot_fish(pop_df, par_df, first_step, last_step, res):
    max_pop = pop_df[["Gen", "Pop"]].groupby("Gen", as_index=False).sum().max()

    parents = par_df["ParentId"].unique()
    children = par_df["ChildId"].unique()
    root_list = np.setdiff1d(parents, children)
    if len(root_list) != 1:
        raise Exception("Failed to determine root. There must be exactly one node with no parent in the parent tree.")
    root_id = root_list[0]
    ids = np.concatenate([[root_id], children])
    tree = build_tree(par_df, ids, root_id)

    # Create Colors
    cols = cm.rainbow(np.linspace(0, 1, len(ids))).tolist()
    random.shuffle(cols)
    col_map = {ids[i]: np.array(list(map(lambda x: int(x * 255), cols[i]))) for i in range(len(ids))}
    col_map[-1] = np.ones(4)  # Root is white

    img_pixels = get_image_pixels(tree, ids, root_id, pop_df, col_map, max_pop, res[1], first_step, last_step, True)
    image = pixels_to_img(img_pixels, res)
    return image


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a Fish (Muller) plot for the given evolutionary tree.')
    parser.add_argument("populations", type=str, help="A CSV file with the header \"Id,Gen,Pop\".")
    parser.add_argument("parent_tree", type=str, help="A CSV file with the header \"ParentId,ChildId\".")
    parser.add_argument("output", type=str, help="Output image filepath. The format must support alpha channels.")
    parser.add_argument("-F", "--first", dest="first_gen", type=int, help="The generation to start plotting from.")
    parser.add_argument("-L", "--last", dest="last_gen", type=int, help="The generation to end the plotting at.")
    parser.add_argument("-W", "--width", dest="width", type=int, help="Output image width", default=1280)
    parser.add_argument("-H", "--height", dest="height", type=int, help="Output image height", default=720)
    args = parser.parse_args()

    resolution = [args.width, args.height]

    populations_df = pd.read_csv(args.populations)
    parent_df = pd.read_csv(args.parent_tree)

    min_gen = populations_df["Gen"].min()
    first_gen = max(min_gen, args.first_gen) if args.first_gen else min_gen

    max_gen = populations_df["Gen"].max()
    last_gen = min(args.last_gen, max_gen) if args.last_gen else max_gen

    if first_gen >= last_gen:
        err = f"First generation {first_gen} must be before the last generation {last_gen}."
        raise Exception(err)

    gen_count = max_gen - min_gen
    populations_df["Step"] = populations_df["Gen"] * resolution[0] / gen_count
    first_step = int(populations_df["Step"].min())
    last_step = int(populations_df["Step"].max())

    img = plot_fish(populations_df, parent_df, first_step, last_step, resolution)
    img.save(args.output)
