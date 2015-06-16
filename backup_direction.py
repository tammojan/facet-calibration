from __future__ import print_function
import numpy as np
import pyrap.tables as pt
import os
import sys
from multiprocessing import Pool


def backup_previous_direction(ms, source, backup=None, column="CORRECTED_DATA"):
    """
    Create a backup of CORRECTED_DATA after running a direction. This data 
    could be entered into SUBTRACTED_DATA_ALL to get the previous state.
    """
    print(ms)
    print(source)
    # Open the MS
    if not os.path.exists(ms):
        raise Exception('MS file not found')
    ms_table = pt.table(ms)
    # Read the column
    if column not in ms_table.colnames():
        raise Exception('Column name not found in the MS')
    data = ms_table.getcol(column)
    # Save the data
    np.save("{}/direction_{}.npy".format(ms, source), data)

    
def restore_direction(ms, backup, column="SUBTRACTED_DATA_ALL"):
    """
    Load a backup of CORRECTED_DATA after running a direction into 
    SUBTRACTED_DATA_ALL.
    """
    print(ms)
    print(backup)
    # Read the data
    data = np.load(backup)
    # Open the MS
    if not os.path.exists(ms):
        raise Exception('MS file not found')
    ms_table = pt.table(ms, readonly=False)
    # Write the column
    if column not in ms_table.colnames():
        raise Exception('Column name not found in the MS')
    ms_table.putcol(column, data)


def backup_previous_direction_p(mslist, source, backup=None, 
                                   column="CORRECTED_DATA", parallel=0):
    if parallel != 0:
        if parallel > 0:
            pool = Pool(parallel)
        else:
            pool = Pool() # Use all the procesors available
        for ms in mslist:
            pool.apply_async(backup_previous_direction, 
                             args=(ms, source,),
                             kwargs={"backup": backup, 
                                     "column": column}
                             )
        pool.close()
        pool.join()
    else:
        for ms in mslist:
            backup_previous_direction(ms, source, backup=backup, 
                                      column=column)



if __name__ == "__main__":
    config = {}
    execfile(sys.argv[1])

    mslist = ["{name:s}_SB{b1:03d}-{b2:03d}.{res:s}.ms".format(name=NAME, res=RES, b1=b, b2=b+9) for b in BANDS]
    
    backup_previous_direction_p(mslist, do_sources[-1], parallel=-1)

