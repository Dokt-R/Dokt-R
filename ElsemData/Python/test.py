import pandas as pd

# pd.set_option('max_rows', 5)

column_names = ['DOY', 'hhmm', 'ssms', 'ch1',
                'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'Battery']

df = pd.read_fwf('../RawData/A3060000.DAT', delimit_whitespace=True,
                 names=column_names)

for col in list(df.columns):
    df[col] = df[col].str[3:]

df['DOY'] = pd.to_numeric(df['DOY'], downcast='integer')
df['Hour'] = df['hhmm'].str[:2]
df['Minute'] = df['hhmm'].str[2:4]
df['Second'] = df['ssms'].str[:2]
df = df.drop(columns=['hhmm', 'ssms', 'Battery'])

df[['ch1', 'ch2', 'ch3', 'ch3', 'ch4', 'ch5', 'ch6']] = df[[
    'ch1', 'ch2', 'ch3', 'ch3', 'ch4', 'ch5', 'ch6']].astype('float32')


print(df)
