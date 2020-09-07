import snakemake
import datetime
import psutil

cores = max(1, psutil.cpu_count() - 1)
stats_file = 'stat_file.txt'

start_time = datetime.datetime.now()

# Please mind dryrun = True!
if not snakemake.snakemake(
        'Snakefile',
        cores=cores,
        printshellcmds=True,
        dryrun=True,
        targets=['all'],
        stats=stats_file):
    raise ValueError("Pipeline failed see Snakemake error message for details")

end_time = datetime.datetime.now()
print("--- Pipeline running time: %s ---" % (str(end_time - start_time)))
