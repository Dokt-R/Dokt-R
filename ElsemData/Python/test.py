import pandas as pd

# pd.set_option('max_rows', 5)

column_names = ['DOY', 'hhmm', 'ssms', 'ch1',
                'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'Battery']

data = pd.read_fwf('../RawData/A3060000.DAT', delimit_whitespace=True,
                   names=column_names)
data = data.drop(columns='Battery')
print(data)
