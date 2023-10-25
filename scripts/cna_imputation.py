import numpy as np
import pandas as pd


def fill_nans(dataframe, chr_lengths):
    grouped = dataframe.groupby(["sampleID", "chrom"]).agg(
        {"start": "min", "end": "max"}
    )
    grouped = grouped.rename(columns={"start": "min_start", "end": "max_end"})
    grouped = grouped.reset_index()
    missing_ranges = []
    for row in grouped.iterrows():
        if row[1].min_start > 1:
            missing_ranges.append(
                [row[1].sampleID, row[1].chrom, 1, row[1].min_start - 1, np.nan, np.nan]
            )
        if row[1].max_end < chr_lengths[f"{row[1].chrom}"] - 1:
            missing_ranges.append(
                [
                    row[1].sampleID,
                    row[1].chrom,
                    row[1].max_end + 1,
                    chr_lengths[f"{row[1].chrom}"],
                    np.nan,
                    np.nan,
                ]
            )
    missing_ends = pd.DataFrame(missing_ranges, columns=dataframe.columns)
    dataframe = pd.concat([dataframe, missing_ends])
    dataframe = dataframe.sort_values(by=["sampleID", "chrom", "start"])
    dataframe.reset_index(inplace=True, drop=True)
    return dataframe


# Makes sure that the columns are of the correct type
def rename_columns(dataframe):
    if dataframe.columns.size != 6:
        raise ValueError(
            "Dataframe must have 6 columns: ",
            "sampleID",
            "chrom",
            "start",
            "end",
            "major_cn",
            "minor_cn",
        )
    dataframe.columns = ["sampleID", "chrom", "start", "end", "major_cn", "minor_cn"]
    # check if any value in column chrom starts with "chr"
    if not dataframe.chrom.str.startswith("chr").any():
        dataframe.chrom = "chr" + dataframe.chrom.astype(str)
    dataframe.sort_values(by=["sampleID", "chrom", "start"], inplace=True)
    dataframe.reset_index(inplace=True, drop=True)
    return dataframe


def are_mergeable(first, second):
    return (
        first.sampleID == second.sampleID
        and first.chrom == second.chrom
        and first.end + 1 == second.start
        and (
            first.minor_cn == second.minor_cn
            or np.isnan(first.minor_cn)
            and np.isnan(second.minor_cn)
        )
        and (
            first.major_cn == second.major_cn
            or np.isnan(first.major_cn)
            and np.isnan(second.major_cn)
        )
    )


def merge_neighbours(sp_df):
    idx_to_remove = []
    to_merge = []
    i = 0
    while i < len(sp_df) - 1:
        first = sp_df.iloc[i]
        second = sp_df.iloc[i + 1]
        if are_mergeable(first, second):
            idx_to_remove.append(i)
            idx_to_remove.append(i + 1)
            # make a copy of first
            merged = first.copy()
            merged.end = second.end
            to_merge.append(merged)
            i += 1
        i += 1

    merged_entries = pd.DataFrame(to_merge)
    # remove from sp_df where idx_to_remove is in the index
    sp_df = sp_df.drop(idx_to_remove)
    # concat the new_table to sp_df
    sp_df = pd.concat([sp_df, merged_entries])
    # sort sp_df by sampleID, chr, start
    sp_df = sp_df.sort_values(by=["sampleID", "chrom", "start"])
    sp_df.reset_index(inplace=True, drop=True)
    return sp_df


def must_merge(sp_df):
    for i in range(0, len(sp_df) - 1):
        first = sp_df.iloc[i]
        second = sp_df.iloc[i + 1]
        if (
            first.sampleID == second.sampleID
            and first.chrom == second.chrom
            and first.end + 1 == second.start
            and first.minor_cn == second.minor_cn
            and first.major_cn == second.major_cn
        ):
            return True
    return False


def merge_loop(data_frame):
    print(f"Entries: {data_frame.shape[0]}")
    while must_merge(data_frame):
        print("merging")
        data_frame = merge_neighbours(data_frame)
        print(f"Entries: {data_frame.shape[0]}. Checking again...")
    print("Finished")
    return data_frame


def create_imputed_entries(dataframe):
    new_entries = []
    for i in range(len(dataframe)):
        if np.isnan(dataframe.loc[i, "major_cn"]):
            if i != 0:
                prev_maj, prev_min = (
                    dataframe.loc[i - 1, "major_cn"],
                    dataframe.loc[i - 1, "minor_cn"],
                )
            if i != len(dataframe) - 1:
                next_maj, next_min = (
                    dataframe.loc[i + 1, "major_cn"],
                    dataframe.loc[i + 1, "minor_cn"],
                )
            if (
                i == 0
                or dataframe.loc[i, "sampleID"] != dataframe.loc[i - 1, "sampleID"]
                or dataframe.loc[i, "chrom"] != dataframe.loc[i - 1, "chrom"]
            ):
                new_entries.append(
                    [
                        dataframe.loc[i, "sampleID"],
                        dataframe.loc[i, "chrom"],
                        dataframe.loc[i, "start"],
                        dataframe.loc[i, "end"],
                        next_maj,
                        next_min,
                    ]
                )
            elif (
                i == len(dataframe) - 1
                or dataframe.loc[i, "sampleID"] != dataframe.loc[i + 1, "sampleID"]
                or dataframe.loc[i, "chrom"] != dataframe.loc[i + 1, "chrom"]
            ):
                new_entries.append(
                    [
                        dataframe.loc[i, "sampleID"],
                        dataframe.loc[i, "chrom"],
                        dataframe.loc[i, "start"],
                        dataframe.loc[i, "end"],
                        prev_maj,
                        prev_min,
                    ]
                )
            else:
                start = dataframe.loc[i, "start"]
                end = dataframe.loc[i, "end"]
                midpoint = start + (end - start) // 2
                new_entries.append(
                    [
                        dataframe.loc[i, "sampleID"],
                        dataframe.loc[i, "chrom"],
                        dataframe.loc[i, "start"],
                        midpoint,
                        prev_maj,
                        prev_min,
                    ]
                )
                new_entries.append(
                    [
                        dataframe.loc[i, "sampleID"],
                        dataframe.loc[i, "chrom"],
                        midpoint + 1,
                        dataframe.loc[i, "end"],
                        next_maj,
                        next_min,
                    ]
                )
    imputataion_df = pd.DataFrame(new_entries, columns=dataframe.columns)
    idx_to_remove = dataframe[dataframe.major_cn.isnull()].index

    
    print(f"New entries: {imputataion_df.shape[0]}")
    print(f"Removed entries: {len(idx_to_remove)}")
    # remove from sp_df where idx_to_remove is in the index
    dataframe = dataframe.drop(idx_to_remove)   
    # concat the new_table to sp_df
    dataframe = pd.concat([dataframe, imputataion_df])
    # sort sp_df by sampleID, chr, start
    dataframe = dataframe.sort_values(by=["sampleID", "chrom", "start"])
    dataframe.reset_index(inplace=True, drop=True)
    return dataframe


def calculate_coverage(dataframe, xx_length, xy_length):
    # Select the rows where major_cn is NaN
    nan_rows = dataframe[~dataframe["major_cn"].isna()]

    # Compute the differences between end and start
    diffs = nan_rows["end"] - nan_rows["start"]

    # Group the differences by sampleID and compute the sum for each group
    sums = diffs.groupby(nan_rows["sampleID"]).sum()

    # Print the resulting Series
    missing_df = pd.DataFrame(sums, columns=["bases"])

    is_male = (
        dataframe.groupby("sampleID")
        .apply(lambda x: (x["chrom"] == "chrY").any())
        .rename("has_Y")
    )

    # merge missing_df with is_male
    missing_df = pd.merge(missing_df, is_male, left_index=True, right_index=True)

    missing_df["frac"] = np.where(
        missing_df["has_Y"],
        missing_df["bases"] / xy_length,
        missing_df["bases"] / xx_length,
    )

    missing_df["chroms"] = dataframe.groupby("sampleID")["chrom"].nunique().tolist()
    return missing_df
