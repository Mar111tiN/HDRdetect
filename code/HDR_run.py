import pandas as pd
from HDR_core import get_HDR_multi, get_filter_hdr_multi
from script_utils_HDR import show_output
from HDR_hotspot import bam2hotspot_multi, pileup2hotspot


HDR_cols = [
    "HDRcand",
    "HDRcount",
    "HDRinfo",
]
base_cols = ["Chr", "Start", "End", "Ref", "Alt"]


def HDR_master(
    mut_file, bam_file, chrom, threads, HDR_config, pileup_file=""
):
    """
    collects the imput and
    """

    
    HDR_config["multi"] = threads > 1
    # Get the mutation file for masterHDR
    show_output(f"Importing {mut_file} for HDR detection", time=False)
    mut_df = pd.read_csv(mut_file, sep="\t").loc[:, base_cols]
    # check for empty file
    if len(mut_df.index) == 0:
        # return an empty df
        show_output(
            f"Mutation file {mut_file} is empty! Writing empty file", color="warning"
        )
        return pd.DataFrame(columns=base_cols)

    # make Chr column categorical for sorting .. and sort
    chrom_list = [f"chr{i}" for i in range(23)] + ["chrX", "chrY"]
    mut_df["Chr"] = pd.Categorical(mut_df["Chr"], chrom_list)
    mut_df = mut_df.sort_values(["Chr", "Start"])

    # output helper
    if HDR_config["multi"]:
        cores = f" using {threads} cores"

    if chrom:
        # output
        show_output(
            f"Starting HDR analysis of {mut_file} on {chrom}{cores}. [MIN_SIM={HDR_config['MINSIM']}, PAD={HDR_config['PAD']}]"
        )
        HDR_df = HDR_run(mut_df, bam_file, chrom, threads, HDR_config, pileup_file)
    else:
        # if no chroms are given
        chroms = mut_df["Chr"].unique()
        HDR_dfs = []
        for chrom in chroms:
            # output
            show_output(
                f"Starting HDR analysis of {mut_file} on {chrom}{cores}. [MIN_SIM={HDR_config['MINSIM']}, PAD={HDR_config['PAD']}]"
            )
            hdr_df = HDR_run(mut_df, bam_file, chrom, threads, HDR_config, pileup_file)
            HDR_dfs.append(hdr_df)

        HDR_df = pd.concat(HDR_dfs).sort_values("Chr")

    return HDR_df


def HDR_run(mut_df, bam_file, chrom, threads, HDR_config, pileup_file):
    """
    does all the work and gets directly called from snakemake/CLI wrapper
    """

    # restrict mut_df to chrom
    mut_df = mut_df.query("Chr == @chrom")
    if mut_df.empty:
        return pd.DataFrame(columns=base_cols + HDR_cols)
    # ########### FIND HOTSPOTS ###########################
    if pileup_file:
        source = pileup_file
        hotspot_df = pileup2hotspot(
            mut_df, pileup_file, chrom=chrom, HDR_config=HDR_config
        )
    else:
        source = bam_file
        # hotspot_df = bam2hotspot(bam_file, chrom=chrom, HDR_config=HDR_config, mut_df=mut_df)
        hotspot_df = bam2hotspot_multi(
            bam_file, chrom=chrom, HDR_config=HDR_config, mut_df=mut_df, threads=threads
        )
    hotspot_len = len(hotspot_df.index)

    # check empty
    if hotspot_len == 0:
        show_output(f"Found no putative HDR lanes in {source}.")
        for col in HDR_cols[:2]:
            mut_df.loc[:, col] = 0
            mut_df.loc[:, 'HDRinfo'] = "no HDR in vincinity"
        return mut_df.loc[:, base_cols + HDR_cols]
        show_output(f"No putative HDR lanes in {source}.")
    show_output(f"Detected {hotspot_len} putative HDR lanes in {source}.")

    # ########### FILTER HOTSPOTS ###########################
    filter_HDR = get_filter_hdr_multi(mut_df, hotspot_df, threads, HDR_config)
    filter_HDR_len = len(filter_HDR.index)

    if filter_HDR_len == 0:
        show_output(f"Found no HDR-rich mutations in {source}", multi=False)
        for col in HDR_cols[:2]:
            mut_df.loc[:, col] = 0
            mut_df.loc[:, 'HDRinfo'] = "no HDR-rich mutations in vincinity"
        return mut_df.loc[:, base_cols + HDR_cols]
    show_output(f"Found a total of {filter_HDR_len} HDR-rich mutations", multi=False)

    HDR_df = get_HDR_multi(
        filter_HDR,
        bam_file=bam_file,
        hotspot_df=hotspot_df,
        threads=threads,
        chrom=chrom,
        HDR_config=HDR_config,
    )

    # merge the HDR output into mut_df
    HDR_df = mut_df.merge(HDR_df, on=base_cols, how="outer").sort_values(base_cols[:2])
    for col in HDR_cols[:2]:
        HDR_df[col] = HDR_df[col].fillna(0).astype(int)
    HDR_df["HDRinfo"] = HDR_df["HDRinfo"].fillna("no HDR in vincinity")

    return HDR_df
