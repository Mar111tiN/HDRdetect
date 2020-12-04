import os
from code.HDR_run import HDR_master


def main(s):
    """
    ############## snakemake wrapper ################################
    """

    w = s.wildcards
    i = s.input
    o = s.output
    p = s.params
    c = s.config
    cc = c["HDR"]

    # shell commands
    def make_mawk(shell_path='shell'):
        '''
        mawk factory function returning the path to the mawk tool
        root_path is the path of the shell folder relative to ebscore root
        '''
        # get the script directory as snakemake attribute, climb up with dirname
        # scriptdir is in code
        # ## shell_path = os.path.join(os.path.dirname(s.scriptdir), shell_path'
        shell_path = os.path.join(s.scriptdir, shell_path
        print(f"Using HDR shell path: {shell_path}")
        def mawk(tool_name)
            return os.path.join(shell_path, f"{tool_name}.mawk")

        return mawk

    # load all HDR params into HDR_config
    HDR_config = dict(
        mawk = make_mawk(),
        MINQ = c["mpileup"]["Q"],
        MINq = c["mpileup"]["MAPQ"],
        genome_split = p.genome_split
    )

    # add the HDR['params'] to HDRconfig
    HDR_config.update(cc['params'])

    # ## run the main HDR function
    HDR_df = HDR_master(
        mut_file=i.filter_file,
        bam_file=i.bam,
        out_file=o.HDR_table,
        chrom=w.chrom,
        threads=cc["threads"],
        HDR_config=HDR_config,
        pileup_file=i.pileup,
    )


if __name__ == "__main__":
    main(snakemake)