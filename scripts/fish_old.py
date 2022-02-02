#!/usr/bin/env python

import argparse
import logging
import random
import sys
from datetime import timedelta
from timeit import default_timer as timer

import numpy as np
import pandas as pd
from matplotlib import cm
from PIL import Image

log = logging.getLogger("Fish.Py")
logging.basicConfig()
log.setLevel(logging.INFO)
eps = sys.float_info.epsilon


# Taken from https://stackoverflow.com/questions/3160699/python-progress-bar
def progressbar(it, prefix="", size=60):
    count = len(it)

    def show(j):
        x = int(size * j / count)
        print("%s[%s%s] %i/%i" % (prefix, "#" * x, "." * (size - x), j, count), end='\r', flush=True)

    show(0)
    for i, item in enumerate(it):
        yield item
        show(i + 1)
    print("\n")


# A dict-based tree where for each parent there is a list of children
def build_tree(par_df, ids, root_id):
    tree = {cid: list(children) for cid in ids if len(children := par_df.loc[par_df["ParentId"] == cid, "ChildId"]) > 0}
    if -1 not in tree:
        tree[-1] = [root_id]
    return tree


# Converts populations to pixels - a color of pixel is a proportional mix of the contributing popluations
def line_from_data(px_row_orig, width):
    result = []
    pop_total = sum(row_it[0] for row_it in px_row_orig)
    # Calculate the fraction of population
    px_row = [[px_row_orig[i][0] * width / pop_total, px_row_orig[i][1]] for i in range(len(px_row_orig))]
    contribution = 1.0
    index = 0
    new_px = np.zeros(4)
    while index < len(px_row):
        if px_row[index][0] >= eps:
            factor = contribution if px_row[index][0] >= contribution else px_row[index][0]
            px_row[index][0] -= factor
            contribution -= factor
            new_px += px_row[index][1] * factor
            if contribution < eps:
                contribution = 1.0
                result.append(new_px)
                new_px = np.zeros(4)
        else:
            index += 1
    if len(result) < width:
        result.append(new_px + 127*contribution*np.ones(4))  # Add a bit of gray to the last pixel if needed
    return result


# Obtains the population by linear interpolation between previous and next know step for the given clone_id
def get_count_at_step(pops_df, first_step, last_step, step, clone_id):
    cp_df = pops_df[pops_df["Id"] == clone_id]
    match = cp_df[cp_df["Step"] == step]
    if not match.empty:
        return float(cp_df["Pop"][match.index])

    l_neigh_i = cp_df[cp_df["Step"] < step]["Step"]
    u_neigh_i = cp_df[cp_df["Step"] > step]["Step"]

    if l_neigh_i.size > 0:
        ind_min = l_neigh_i.idxmax()
        prev_step, prev_pop = cp_df.loc[ind_min, ["Step", "Pop"]].tolist()
    else:
        prev_step, prev_pop = first_step, 0

    if u_neigh_i.size > 0:
        ind_max = u_neigh_i.idxmin()
        next_step, next_pop = cp_df.loc[ind_max, ["Step", "Pop"]].tolist()
    else:
        next_step, next_pop = last_step, 0

    if prev_pop + next_pop == 0:
        return 0

    pop = ((step - prev_step) * next_pop + (next_step - step) * prev_pop) / (next_step - prev_step)
    log.debug(f"Id: {clone_id}, Step: {step}, Pop: {pop}")
    return pop


# Population at step for each clone ind ids, if normalized the remainder is given to -1
def get_pop_dict(pops_df, ids, fg, lg, step, max_pop, norm):
    pop_dict = {clone_id: count for clone_id in ids if
                (count := get_count_at_step(pops_df, fg, lg, step, clone_id)) > 0}
    if not norm:
        pop_count = sum(pop_dict.values())
        remainder = max_pop - pop_count
        pop_dict[-1] = remainder
    return pop_dict


# Construct the list of the form [(id, pop_size),...] by inserting children's populations in the middle of parent's
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


def get_image_pixels(tree, ids, root_id, pop_df, col_map, max_pop, breath, first_step, last_step, norm):
    img_pixels = []
    for step in progressbar(range(first_step, last_step + 1), "Step: "):
        pop_dict = get_pop_dict(pop_df, ids, first_step, last_step, step, max_pop, norm)
        if len(pop_dict) > 0:
            id_rows = rec_descend(tree, pop_dict, root_id if norm else -1)
            pixel_row = [(id_row[0], col_map[id_row[1]]) for id_row in id_rows]
        else:
            pixel_row = [(-1, np.ones(4) * 127)]  # use gray if there's no population
        line_px = line_from_data(pixel_row, breath)
        img_pixels.append(line_px)
    return img_pixels


def plot_fish(pop_df, par_df, first_step, last_step, res, norm):
    max_pop = int(pop_df[["Step", "Pop"]].groupby("Step", as_index=False).sum().max()["Pop"])

    # Build the internal search tree
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

    # Convert
    img_pixels = get_image_pixels(tree, ids, root_id, pop_df, col_map, max_pop, res[1], first_step, last_step, norm)
    image = pixels_to_img(img_pixels, res)
    return image


if __name__ == '__main__':
    randint = random.randint(0, sys.maxsize)
    parser = argparse.ArgumentParser(description='Create a Fish (Muller) plot for the given evolutionary tree.')
    parser.add_argument("populations", type=str, help="A CSV file with the header \"Id,Step,Pop\".")
    parser.add_argument("parent_tree", type=str, help="A CSV file with the header \"ParentId,ChildId\".")
    parser.add_argument("output", type=str, help="Output image filepath. The format must support alpha channels.")
    parser.add_argument("-F", "--first", dest="first_step", type=int, help="The step to start plotting from.")
    parser.add_argument("-L", "--last", dest="last_step", type=int, help="The step to end the plotting at.")
    parser.add_argument("-W", "--width", dest="width", type=int, help="Output image width", default=1280)
    parser.add_argument("-H", "--height", dest="height", type=int, help="Output image height", default=720)
    parser.add_argument("-S", "--seed", dest="seed", type=int, help="Seed for colors", default=randint)
    parser.add_argument('-r', dest="raw", default=False,
                        help='Plot the individual steps rather than interpolated values', action='store_true')
    parser.add_argument('-a', dest="absolute", default=False,
                        help='Plot the populations in absolute numbers rather than normalized', action='store_true')

    args = parser.parse_args()

    resolution = [args.width, args.height]
    random.seed(args.seed)
    log.info(f"Creating {args.width}*{args.height} Fish plot file {args.output} with the seed {args.seed}")

    populations_df = pd.read_csv(args.populations)
    parent_df = pd.read_csv(args.parent_tree)

    min_step = populations_df["Step"].min()
    fs = max(min_step, args.first_step) if args.first_step else min_step

    max_step = populations_df["Step"].max()
    ls = min(args.last_step, max_step) if args.last_step else max_step

    if fs >= ls:
        err = f"First step {fs} must be before the last step {ls}."
        raise Exception(err)

    populations_df = populations_df[(populations_df["Step"] >= fs) & (populations_df["Step"] <= ls)]

    if not args.raw:
        step_count = max_step - min_step
        norm_step_count = resolution[0] - 1
        populations_df["Step"] = populations_df["Step"] * norm_step_count / step_count
        fs = int(populations_df["Step"].min())
        ls = fs + norm_step_count

    normalize = not args.absolute
    log.info("Start timer")
    start = timer()
    img = plot_fish(populations_df, parent_df, fs, ls, resolution, normalize)
    end = timer()
    log.info(timedelta(seconds=end-start))
    img.save(args.output)
