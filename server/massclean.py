import logging
import timeit
from os import walk, mkdir
from os.path import exists, join, splitext
import pandas as pd

# Logger Setup

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s:%(levelname)s -> %(message)s')

file_handler = logging.FileHandler('clean.log', mode='w')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)

# Start timer and declare constant variables
starttime = timeit.default_timer()

PATH = r".\dat"
COLUMNS = ['ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6']

for subdir, dirs, files in walk(PATH):
    for file in files:
        filePath = join(subdir, file)
        csvName = splitext(file)[0]+'.csv'

        newDir = subdir.replace('dat', 'csv')
        if not exists(newDir):
            mkdir(newDir)
        csvPath = join(newDir, csvName)

        try:
            codename = file[0]
            year = int(file[1:5])
            doy = int(file[5:8])
        except ValueError as e:
            logger.error('%s\nLocation: %s', e, filePath)
        except Exception:
            logger.exception('Something unexpected happened!')
        else:
            date = pd.to_datetime(year * 1000 + doy, format='%Y%j')
            daterange = pd.date_range(date, periods=86400, freq='s')
            df = pd.read_csv(
                filePath, usecols=[4, 5, 6, 7, 8, 9], names=COLUMNS)
            df.insert(loc=0, column='moment', value=daterange)
            df.insert(loc=0, column='codename', value=codename)

            df.to_csv(csvPath, index=False)

            logger.info('%s finished loading', filePath)

        break

endtime = timeit.default_timer()
print(endtime-starttime, "s")
