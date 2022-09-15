# Elsem-Net Data Management Project

<div style="text-align: justify">
<a href="http://elsem-net.uniwa.gr/">Elsem-Net</a> is a network of telemetric
stations spanning all over Greece for the monitoring of fracture induced
electromagnetic emissions. This network produces daily files for each station
which in turn are sent to a Network Attached Storage.
</div>

<small>_*This project along with its documentation is currently a work in progress_</small>

- [X] Folders
- [X] Current File Structure
- [ ] Legacy File Structure
- [ ] Analysis Files
- [ ] Automation

## The Problem

<div style="text-align: justify">
When starting my diploma thesis I never imagined the scope of the project. It
was supposed to be a simple website utilizing a NAS for its database. What I
didn't know was that after finishing the thesis I would be leaving behind
something incomplete. This network produces 86400 measurements daily across 6
channels per station, 10 stations total. On top of that there are about 30 years
of legacy data that needed to be cleaned and stored in a way that makes it easy
to access, query and analyze.</div>

## General Information

### Folder Structure

<div style="text-align: justify">
The folder structure of <i>most</i> of the data is neatly organized. First by
year and then by station.
</div>

```
Data
└───2015
|   └───A
|   |   └───A2015001.DAT
|   |       A2015002.DAT
|   |       ...
|   └───J
|   |   └───...
|   |       J2015364.DAT
|   |       J2015365.DAT
|   ...
└───2016
|   └───...
...
```

### File Naming

<div style="text-align: justify">
<strong><em>Daily</em></strong> files follow the naming convention bellow without
exception.
</div>

|Information        | Filename       |
|-------------------|----------------|
|Station Codename   |***A*** 2021036.DAT|
|Year               |A ***2021*** 036.DAT|
|Day of Year        |A2021 ***036*** .DAT|

### Daily Data

<div style="text-align: justify">
All legacy data, which is the data that this project currently focuses on comes
in a form different than hourly data. It is an already merged file from all the
hourly files that came from the dataloggers in a zipped format and some cleaning
has been done. However it still has some information that we do not need to store
along with a hard to process time series format.
</div>

```
102,0001.,0000.,00.00,04.93,00.00,212.1,23.43,573.3,660.9,13.64
102,0001.,0000.,01.00,06.17,00.00,207.1,23.43,567.2,662.1,13.64
102,0001.,0000.,02.00,04.93,00.00,202.2,22.19,569.6,676.9,13.64
102,0001.,0000.,03.00,06.17,00.00,217.0,23.43,572.1,679.4,13.64
102,0001.,0000.,04.00,06.17,00.00,207.1,22.19,567.2,678.2,13.64
...
```

### Hourly Data

<div style="text-align: justify">
The data gathered by the dataloggers are collected in this form as hourly data
and I will get to this when the project is up to date. Legacy data are a priority
and are already converted to the daily form mentioned above.<br>Hourly data
looks like this:
</div>

```
...
02+0306.  03+0000.  04+58.00  05+00.00  06+00.00  07+197.3  08+25.89  09+562.2  10+731.2  11+13.64
02+0306.  03+0000.  04+59.00  05+00.00  06+00.00  07+180.0  08+24.66  09+558.5  10+710.2  11+13.64
02+0306.  03+0001.  04+00.00  05+00.00  06+00.00  07+175.1  08+25.89  09+574.6  10+731.2  11+13.64
02+0306.  03+0001.  04+01.00  05+00.00  06+00.00  07+183.7  08+24.66  09+556.1  10+721.3  11+13.64
02+0306.  03+0001.  04+02.00  05+01.00  06+00.00  07+203.4  08+24.66  09+549.9  10+707.7  11+13.64
...
```
<small>_*Needs example with a better ch6 value_</small>

Having some basic information about the stations (86400 rows, 6 channels, Battery Measurements)
we can see that the data follows the pattern below:

```
DayOfYear, hhmm, ss.ms, ch1, ch2, ch3, ch4, ch5, ch6, Battery
```

Therefore:
```
column_names = ['DOY', 'hhmm', 'ssms', 'ch1',
                'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'Battery']
```
Which can be read with pandas as a csv using whitespace delimiter:

```
pd.read_csv('../RawData/A3060000.DAT', delim_whitespace=True, names = column_names)
```

Alternatively, since technically it is not a .csv but a fixed width .DAT file:

```
pd.read_fwf('../RawData/A3060000.DAT')
```
<small>_*Needs more testing_</small>

## Data Cleaning

<div style="text-align: justify">
It is obvious that our data has information in front of every row that makes it
hard to read data and is not useful for analysis. This has a similar pattern
of leading 3 characters, and they can be easily removed.
</div>

```
for col in list(df.columns):
    df[col] = df[col].str[3:]
```

### DOY Column

<div style="text-align: justify">
Days of year also follows a pattern, leading zeroes, the day of year and a decimal.
Instead of using a regular expression however we can simply convert the string to
an integer, since no days will ever be float and save a step by converting it to
the correct data type.
</div>

```
df['DOY'] = pd.to_numeric(df['DOY'], downcast='integer')
```

### Hours, minutes, seconds

### Channels

Every channel is a float data type and is converted using the same logic as DOY.

```
df[['ch1', 'ch2', 'ch3', 'ch3', 'ch4', 'ch5', 'ch6']] = df[[
    'ch1', 'ch2', 'ch3', 'ch3', 'ch4', 'ch5', 'ch6']].astype('float32')
```
<div style="text-align: justify">
<small>_The astype() method could be used from the start, adding int type for DOY
            using only 1 line of code. Two different techniques were used mainly to 
            explore the options that one has when working with data._</small>
</div>

### Battery Column

Since battery values will never be used we can drop the column using:

```
data = data.drop(columns='Battery')
```