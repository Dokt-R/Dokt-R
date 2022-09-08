### Elsem-Net Database Project

This projext along with its documentation is currently a work in progress

## Problem

Elsem-Net is a network of telemetric stations spanning all over Greece for the monitoring of
freacture induced electromagnetic emissions. This network produces daily files for each station
which in turn is sent to a Network Attached Storage. This is in the form of zipped files
containing several fixed width .DAT hourly files that follow the naming convetion bellow

|Information        | Filename       |
|-------------------|----------------|
|Station Codename   |***A***2021036.DAT|
|Year               |A***2021***036.DAT|
|Day of Year        |A2021***036***.DAT|

## Loading Data
The data gathere by the dataloggers are collected in this form as houtly data:

02+0306.  03+0000.  04+58.00  05+00.00  06+00.00  07+197.3  08+25.89  09+562.2  10+731.2  11+13.64
02+0306.  03+0000.  04+59.00  05+00.00  06+00.00  07+180.0  08+24.66  09+558.5  10+710.2  11+13.64
02+0306.  03+0001.  04+00.00  05+00.00  06+00.00  07+175.1  08+25.89  09+574.6  10+731.2  11+13.64
02+0306.  03+0001.  04+01.00  05+00.00  06+00.00  07+183.7  08+24.66  09+556.1  10+721.3  11+13.64
02+0306.  03+0001.  04+02.00  05+01.00  06+00.00  07+203.4  08+24.66  09+549.9  10+707.7  11+13.64
!!! Needs example with a better ch6 value

Having some basic information about the stations (86400 rows, 6 channels, Battery Measurements)
we can see that the data follows the pattern below:

DayOfYear, hhmm, ss.ms, ch1, ch2, ch3, ch4, ch5, ch6, Battery

Therefore:
column_names = ['DOY', 'hhmm', 'ssms', 'ch1',
                'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'Battery']

Which can be read with pandas as a csv using whitespace delimiter:
pd.read_csv('../RawData/A3060000.DAT', delim_whitespace=True, names = column_names)

Alternatively, since technically it is not a .csv but a fixed width .DAT file:
pd.read_fwf('../RawData/A3060000.DAT')
!!! Needs testing for values of 1000+

## Data Cleaning

It is obvious that our data has information in front of every row that makes it
hard to read data and is not useful for analysis. This has a similar pattern
of leading 3 characters, and they can be easily removed.

for col in list(df.columns):
    df[col] = df[col].str[3:]

### DOY Column

Days of year also follows a pattern, leading zeroes, the day of year and a decimal.
Instead of using a regular expression however we can simply convert the string to
an integer, since no days will ever be float and save a step by converting it to
the correct data type.

df['DOY'] = pd.to_numeric(df['DOY'], downcast='integer')

### Hours, minutes, seconds

### Channels

Every channel is a float data type and is converted using the same logic as DOY.

df[['ch1', 'ch2', 'ch3', 'ch3', 'ch4', 'ch5', 'ch6']] = df[[
    'ch1', 'ch2', 'ch3', 'ch3', 'ch4', 'ch5', 'ch6']].astype('float32')

/Comment:   The astype() method could be used from the start, adding int type for DOY
            using only 1 line of code. Two different techniques were used mainly to 
            explore the option that one has when working with data.

### Battery Column

Since battery values will never be used we can drop the column using:
data = data.drop(columns='Battery')
