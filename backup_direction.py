from __future__ import print_function
import numpy as np
import pyrap.tables as pt
import os
from multiprocessing import Pool


def backup_previous_direction(ms, source, backup=None, column="CORRECTED_DATA"):
    """
    Create a backup of CORRECTED_DATA after running a direction. This data 
    could be entered into SUBTRACTED_DATA_ALL to get the previous state.
    """
    # Open the MS
    if not os.path.exists():
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
    # Read the data
    data = np.load(backup)
    # Open the MS
    if not os.path.exists():
        raise Exception('MS file not found')
    ms_table = pt.table(ms, readonly=False)
    # Write the column
    if column not in ms_table.colnames():
        raise Exception('Column name not found in the MS')
    ms_table.putcol(column, data)

    
if __name__ == "__main__":
    config = {}
    execfile(sys.argv[1])

    mslist = ["{name:s}_SB{b1:03d}-{b2:03d}.{res:s}.ms".format(name=NAME, res=RES, b1=b, b2=b+9) for b in BANDS]
    
    pool = Pool()
    for ms in mslist:
        pool.apply_async(backup_previous_direction, ms, do_sources[-1])
    pool.close()
    pool.wait()
