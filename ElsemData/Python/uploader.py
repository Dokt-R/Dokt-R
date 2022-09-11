''' This script is used after running converter.
It walks through all folders in a given path,
and uploads them to the database specified using
LOAD DATA INFILE. This allows faster uploading than
a regular line by line insert but assumes that .csv
files are already in the correct format.'''

import logging
import os
import mysql.connector

# Logger Settings
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s -> %(message)s')
file_handler = logging.FileHandler('uploader.log', 'w')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

# Connections Settings
conn = mysql.connector.connect(
    host="localhost", user="root", password="", database="elsemdb")

# Walking Path and LOAD INFILE
PATH_CSV = os.environ.get('PATH_CSV')
PATH_SQL = os.environ.get('PATH_SQL')

for dirs, subdir, files in os.walk(PATH_CSV):
    for file in files:
        # SQL runs into errors if \ is used so it must be converted
        filepath = os.path.join(PATH_SQL, dirs, file).replace(os.sep, '/')
        sql = f"LOAD DATA INFILE '{filepath}' IGNORE INTO TABLE test FIELDS TERMINATED BY ','"
        try:
            cursor = conn.cursor(buffered=True)
            cursor.execute(sql)
            conn.commit()
        except Exception as e:
            # As of this point I do not know enough about exceptions
            # I haven't encountered a specific error so a generic
            # try/except is used for logging.
            logger.error("File %s did not upload\n%s", filepath, e)
        else:
            logger.info("%s file finished loading successfully", filepath)
        finally:
            cursor.close()
conn.close()
