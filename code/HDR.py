import os
from HDR_core import run_HDR


def main(s):
    # ############## s ################################
    w = s.wildcards
    config = s.config
    i = s.input
    p = s.params
    PAD = min(config['HDR']['padding'], config['filter_bam']['padding'])
    threads = config['HDR']['threads']
    bam = os.path.split(s.input.filter_bam)[1]
    bam = os.path.join(config['filter_bam']['folder'], bam)

    # get the tumor and normal bam files from the tumor-normal file name
    # remove the normal part of the tumor-normal descriptor
    tumor_bam = bam.replace(
        f"-{w.normal}", '').replace('.done', '.bam').replace('filterbamdone', 'filterbam')
    # remove the tumor part of the tumor-normal descriptor
    normal_bam = bam.replace(
        f"{w.tumor}-", '').replace('.done', '.bam').replace('filterbamdone', 'filterbam')
    # print('tumor:', tumor_bam, 'normal:', normal_bam)

    # ## run the main HDR function
    run_HDR(
        mut_file=i.filter_file,
        tumor_bam=tumor_bam,
        normal_bam=normal_bam,
        chrom=w.chrom,
        filter_pileup=i.pileup,
        out_file=s.output.HDR_table,
        threads=threads,
        MINSIM=p.min_sim,
        PAD=PAD,
        MINQ=p.min_q,
        HDRMINCOUNT=p.min_HDR_count
    )


if __name__ == "__main__":
    main(snakemake)


def run_HDR(mut_file, tumor_bam, normal_bam, chrom, filter_pileup, out_file, threads, MINSIM=.90, PAD=100, MINQ=25, HDRMINCOUNT=1, MINq=20):
    '''
    runs the HDR computation for both tumor and normal and populates the table with the respective columns
    '''

    # run HDR for both tumor and normal
    bam_dict = {'Tumor': tumor_bam, 'Normal': normal_bam}

    # ###### PILEUP ANALYSIS ##############################
    for T_or_N in bam_dict.keys():
        show_output(f'Analysing {T_or_N}')
        HDR_df = HDR_master(bam_dict[T_or_N], HDR_df, chrom=chrom, pileup_file=filter_pileup, threads=threads, _type=T_or_N,
                            HDR_config['MINSIM']=HDR_config['MINSIM'], padding=PAD, min_HDR_count=HDRMINCOUNT, MINQ=MINQ, MINq=MINq)

        HDR_df = HDR_df.rename(columns={
            'HDR': f'{T_or_N}HDRcand',
            'HDRcount': f'{T_or_N}HDRcount',
            'HDRinfo': f'{T_or_N}HDRinfo'
        })

    HDR_len = len(HDR_df.query('TumorHDRcount > 0').index)

    show_output(
        f"Found {HDR_len} possible HDR mutations. Writing to file {out_file}")
    HDR_df.to_csv(out_file, sep='\t', index=False)
