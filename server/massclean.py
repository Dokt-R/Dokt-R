import logging
import timeit
from os import listdir, walk
from os.path import isfile, join, splitext
import pandas as pd

logging.basicConfig(filename="log.txt", level=logging.DEBUG,
                    format="%(asctime)s %(message)s", filemode='w')

starttime = timeit.default_timer()


class FileProperties:
    def __init__(self, file) -> None:
        self.filename = file.split('.')[0]
        self.codename = file[0]
        self.year = int(file[1:5])
        self.doy = int(file[5:8])
        self.date = pd.to_datetime(self.year * 1000 + self.doy, format='%Y%j')
        self.daterange = pd.date_range(self.date, periods=86400, freq='s')


# Script


PATH = ".\\dat"
PATH_CSV = ".\\csv"

filelist = [f for f in listdir(PATH) if isfile(join(PATH, f))]
COLUMNS = ['ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6']

for subdir, dirs, files in walk(PATH):
    for file in files:
        print(join(subdir, file))
        print(file)
        # f = FileProperties(file)

        # df = pd.read_csv(
        #     PATH+file, usecols=[4, 5, 6, 7, 8, 9], names=COLUMNS)

        # df.insert(loc=0, column='moment', value=f.daterange)
        # df.insert(loc=0, column='codename', value=f.codename)

        # df.to_csv(PATH_CSV+filename+'.csv', index=False)
        # print(f.filename, f.year, f.doy, f.date)
        # print(df)
        break

endtime = timeit.default_timer()
print(endtime-starttime, "s")
