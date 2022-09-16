# Elsem-Net Data Management Project

*<a href="http://elsem-net.uniwa.gr/">Elsem-Net</a>* is a network of telemetric
stations spanning all over Greece for the monitoring of fracture induced
electromagnetic emissions. This network produces daily files for each station
which in turn are sent to a Network Attached Storage.

<small>_*This project along with its documentation is currently a work in progress_</small>

- [X] Folders
- [X] Current File Structure
- [ ] Legacy File Structure
- [ ] Analysis Files
- [ ] Automation

## The Problem

When starting my diploma thesis I never imagined the scope of the project. It
was supposed to be a simple website utilizing a NAS for its database. What I
didn't know was that after finishing the thesis I would be leaving behind
something incomplete. This network produces 86400 measurements daily across 6
channels per station, 10 stations total. On top of that there are about 30 years
of legacy data that needed to be cleaned and stored in a way that makes it easy
to access, query and analyze.

## General Information

### Folder Structure

The folder structure of *most* of the data is neatly organized. First by
year and then by station.

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

***Daily*** files follow the naming convention bellow without
exception.

|Information        | Filename       |
|-------------------|----------------|
|Station Codename   |***A*** 2021036.DAT|
|Year               |A ***2021*** 036.DAT|
|Day of Year        |A2021 ***036*** .DAT|

### Daily Data

All legacy data, which is the data that this project currently focuses on comes
in a form different than hourly data. It is an already merged file from all the
hourly files that came from the dataloggers in a zipped format and some cleaning
has been done. However it still has some information that we do not need to store
along with a hard to process time series format.

```
102,0001.,0000.,00.00,04.93,00.00,212.1,23.43,573.3,660.9,13.64
102,0001.,0000.,01.00,06.17,00.00,207.1,23.43,567.2,662.1,13.64
102,0001.,0000.,02.00,04.93,00.00,202.2,22.19,569.6,676.9,13.64
102,0001.,0000.,03.00,06.17,00.00,217.0,23.43,572.1,679.4,13.64
102,0001.,0000.,04.00,06.17,00.00,207.1,22.19,567.2,678.2,13.64
...
```

### Hourly Data

The data gathered by the dataloggers are collected in this form as hourly data
and I will get to this when the project is up to date. Legacy data are a priority
and are already converted to the daily form mentioned above.<br>Hourly data
looks like this:

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

## Legacy Files Data Cleaning

The [converter](https://github.com/Dokt-R/Dokt-R/blob/main/ElsemData/Python/converter.py)
file does all the cleaning. It opens every folder, mirrors the directory into a
./csv/ path and then converts all .DAT files to the final .csv that will be loaded
into the database. Saving a csv instead of working on the files as is was the best
solution, as uploading thousands of rows with pandas .to_sql() method was much 
slower, and a sacrifice in space was a small compromise compared to the time it 
saves. On top of that the .csv format is kept as a backup and used in MATLAB for
analysis.

### Opening Files

Since legacy files are already sorted neatly into folders all we have to do is
iterate through them. This two lines make it a breeze:

```
for subdir, dirs, files in os.walk(PATH):
    for file in files:
```

### Database Columns

Each measurement can be easily defined by *two values*. The *moment* that the 
measurement was made and the *station* that made it. Since files are named after 
the station's *codename* it makes sense to have the first column as a **FOREIGN
KEY** with that *codename*. The second column is the moment in datetime to make it 
easier to plot. The order of this is also important since by the end of the project
there will be billion rows of data we want an ***index*** that is efficient.
[This website](https://use-the-index-luke.com/sql/where-clause/searching-for-ranges/greater-less-between-tuning-sql-access-filter-predicates) explains why **index for
equality should be used first.**

Now there are plenty solutions that can get the right daterange for each file.
We already mentioned that the measurements are taken each second but the formating
was ugly, using multiple columns to define one point in time. Instead of reshaping
everything, a pandas daterange from scratch seemed faster. All we have to use is
the filename's Day Of Year.

```
pd.to_datetime(year * 1000 + doy, format='%Y%j')
daterange = pd.date_range(date, periods=86400, freq='s')
```

The rest of the columns are loaded using ```.read_csv(filePath, usecols=[4, 5, 6, 7, 8, 9])``` method.

## Hourly Files Data Cleaning

It is obvious that hourly data has information in front of every row that makes it
hard to read data and is not useful for analysis. This has a similar pattern
of leading 3 characters, and they can be easily removed.

```
for col in list(df.columns):
    df[col] = df[col].str[3:]
```

### DOY Column

Days of year also follows a pattern, leading zeroes, the day of year and a decimal.
Instead of using a regular expression however we can simply convert the string to
an integer, since no days will ever be float and save a step by converting it to
the correct data type.

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
<small>_The astype() method could be used from the start, adding int type for DOY
            using only 1 line of code. Two different techniques were used mainly to 
            explore the options that one has when working with data._</small>

### Battery Column

Since battery values will never be used we can drop the column using:

```
data = data.drop(columns='Battery')
```