'''This script walks through given directory
and converts all files to .csv clean data files
ready to be imported to MariaDB.'''
import logging
import os
import pandas as pd

# Logger Setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s:%(levelname)s -> %(message)s')
file_handler = logging.FileHandler('clean.log', mode='w')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

# Variables and Converter
PATH = r".\dat"
COLUMNS = ['ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6']

for subdir, dirs, files in os.walk(PATH):
    for file in files:
        filePath = os.path.join(subdir, file)
        csvName = os.path.splitext(file)[0]+'.csv'

        # Ensure that the same folder structure is used
        newDir = subdir.replace('dat', 'csv')
        if not os.path.exists(newDir):
            os.mkdir(newDir)
        csvPath = os.path.join(newDir, csvName)

        try:
            codename = file[0]
            year = int(file[1:5])
            doy = int(file[5:8])
        except ValueError as e:
            logger.error('%s\nLocation: %s', e, filePath)
        except Exception as e:
            logger.exception('Something unexpected happened! \n%s', e)
        else:
            date = pd.to_datetime(year * 1000 + doy, format='%Y%j')
            daterange = pd.date_range(date, periods=86400, freq='s')
            df = pd.read_csv(
                filePath, usecols=[4, 5, 6, 7, 8, 9], names=COLUMNS)
            df.insert(loc=0, column='moment', value=daterange)
            df.insert(loc=0, column='codename', value=codename)

            df.to_csv(csvPath, index=False)

            logger.info('%s finished loading', filePath)
