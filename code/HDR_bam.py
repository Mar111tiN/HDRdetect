import os
import pandas as pd
import sys
from script_utils_HDR import cmd2df


def bam2df2(
    bam_file,
    HDR_config,
    bamtags="",
    tool="samtools",
    mut_row=None,
    region="",
):
    """
    set the region requires 3 threads
    """

    # unwrap the mawk tools using mawk tool unwrapper
    def mawk(tool):
        return os.path.join(HDR_config["mawk_path"], f"{tool}.mawk")

    if region:
        mut_pos = region
    elif isinstance(mut_row, pd.Series):
        chrom = mut_row["Chr"]
        pos = mut_row["Pos"]
        mut_pos = f"{chrom}:{pos}-{pos} "
        print(mut_pos)
    else:
        mut_pos = ""
    cmd = (
        f"{tool} view {bam_file} {mut_pos} | {mawk('bam2csv')} | {mawk('editbam')} {HDR_config['MINq']}"
    )
    # show_command(cmd)
    bam_df = cmd2df(cmd, show=False)
    return bam_df
