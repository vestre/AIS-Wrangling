"""
File for executing whole AIS workflow
"""

# Importing relevant libraries
import numpy as np
import pandas as pd 
import numba
import datetime

from NearCollisionStatistics *

### Reading from file and sorting ###

# Setting internal options for speedup
pd.set_option('compute.use_bottleneck', True)
pd.set_option('compute.use_numexpr', True)

# Setting the path
PATH = "~/DATA/aisdata.csv"

# Setting columns to read
incols = ["date_time_utc","mmsi","lat","lon","sog","cog","true_heading",
            "length","nav_status","RISK_Norwegian_Main_Vessel_Category_ID"]

# Setting datecolumn if we want to convert times directly
datecol = 2

AIS_df = pd.io.parsers.read_csv(PATH
                ,engine="c"
                ,sep=";"
                ,usecols=incols        
                ,compression=None
                ,dtype={"mmsi": np.uint64,             
                        "lat": np.float32, 
                        "lon": np.float32,
                        "sog": np.float32,
                        "cog": np.float32,
                        "true_heading": np.float32,
                        "length": np.float32,
                        "nav_status": np.int32,
                        "RISK_Norwegian_Main_Vessel_Category_ID": np.float32}
                        )

# Renaming fields
AIS_df.rename(columns={"date_time_utc": "Time",
                        "mmsi": "ID",
                        "lon": "LON",
                        "lat": "LAT",
                        "nav_status": "Status",
                        "true_heading": "Heading",
                        "sog": "SOG",
                        "cog": "COG",
                        "length": "Length",
                        "RISK_Norwegian_Main_Vessel_Category_ID": "Category"
                         },inplace=True)

# Converting time column from string to datetime
AIS_df["Time"] = pd.to_datetime(AIS_df.Time,
                            format="%Y-%m-%dT%H:%M:%S.%f",box=False)

# Storing datetimes in the dataframe
AIS_df["Time_datetime"] = AIS_df["Time"]

# Converting datetimes to ints for later conversion to intervals
AIS_df["Time"] = AIS_df.Time.values.astype(np.uint64) / 10**9

# Removing platforms
AIS_df = AIS_df[(AIS_df.ID > 99999999)&(AIS_df.ID <= 999999999)]

# Sorting by time
AIS_df.sort_values(by = "Time",inplace=True)

### Making intervals on the sotred table and filtering for unique vessels ###

# Making a new pointer to the table
AIS_sorted = AIS_df

# Creating intervals of 20 seconds
AIS_sorted["Interval"] = (AIS_sorted.Time - np.min(AIS_sorted.Time)) // 20

# Optional: Creating ntervals of 60 minutes
AIS_sorted["IntervalMins"] = (AIS_sorted.Time - np.min(AIS_sorted.Time)) // 60

# Making Interval iterable (to integer values)
AIS_sorted.Interval = AIS_sorted.Interval.values.astype(np.uint32)
AIS_sorted.IntervalMins = AIS_sorted.IntervalMins.values.astype(np.uint32)

# Dropping duplicates if they have both the same ID and the same interval, 
# except for the first unique instance - going from 18 to 8 M instances
# This equates to finding each unique instance in each interval
AIS_shorter = AIS_sorted.drop_duplicates(["ID","Interval"],keep="first")

# Sorting in order to speed up computations (making LAT and LONs in increasing
# order (this is important for the Time to CPA calculations))
# NB: Is now implemented in the function as well.
AIS_shorter = AIS_shorter.sort_values(by=["Time","LAT","LON"], 
                                        axis=0, inplace=False)

### Calculating time to CPA and CPA distance ###

# Calculating times to CPA and CPA distances for all values in all
# intervals
# NB: This is the most computationally heavy operation - see function code
AIS_outmatrix = time_to_CPA_calculator(AIS_shorter, extrainterval=False)

### Filtering the data by tCPA and CPA ###

# Making pointer to the output matrix
AIS_filtered = AIS_outmatrix

# Filtering by SOG and Time to CPA and Navigation status in order to exclude
# situations which are not in fact near-collision situations, and where the
# vessels are in tugging mode (11, 12 status), as well as excluding vessels
# that are not moving.
SOG_tresh_min = 5
SOG_tresh_max = 50
filter3 = (((AIS_outmatrix["SOG_Vessel_1"]>SOG_tresh_min)
            &(AIS_outmatrix["SOG_Vessel_2"]>SOG_tresh_min))
            &((AIS_outmatrix["SOG_Vessel_1"]<SOG_tresh_max)
            &(AIS_outmatrix["SOG_Vessel_2"]<SOG_tresh_max))
            &((AIS_outmatrix["Status_Vessel_1"] != 11)
            &(AIS_outmatrix["Status_Vessel_1"] != 12))
            &((AIS_outmatrix["Status_Vessel_2"] != 11)
            &(AIS_outmatrix["Status_Vessel_2"] != 12))
            &(AIS_outmatrix["Min_time_to_CPA"] < 1000)
            &(AIS_outmatrix["CPA"] < 1000))

# Applying the filter
AIS_filtered_more = AIS_outmatrix[filter3]

# Sorting the table in order to later drop only high values of min time to
# CPA
AIS_Filtered_more = AIS_filtered_more.sort_values(["Min_time_to_CPA"])

# Finding a unique record representing each near-collision siutation, choosing
# the record with the lowest time to CPA (by the sorting above)
AIS_Unique_Filtered_more = AIS_Filtered_more.drop_duplicates(
                                        ["ID_Vessel_1","ID_Vessel_2"])

# Sorting further by minimum time to CPA (this might be superflous)
AIS_Unique_Filtered_more = AIS_Unique_Filtered_more.sort_values(
                                        ["Min_time_to_CPA"])

### Building a databank ###

# Renaming pointer
AIS_filtered = AIS_Unique_Filtered_more

# Collecting observations by using the filtered list of situations as well
# as the original databank.
DATABANK = observation_collector(AIS_df, AIS_filtered)

# Releasing the memory used by the large AIS data frame and the smaller filter
AIS_df = 0
AIS_sorted = 0
AIS_shorter = 0
AIS_filtered = 0
AIS_outmatrix = 0
AIS_Unique_Filtered_more = 0

### Synchronizing time scales in databank ###

# Creating pointer
AIS_databank = DATABANK

# Synchronizing times
databank_mergetime = observation_synchronizer(AIS_databank)

### Calculating Statistics (CPA, tCPA, CPA_dist) for all observations ###

# Creating pointer
AIS_db = databank_mergetime

# Calculating stats
AIS_with_stats = Stat_calculator(AIS_db)

# Calculating dSOG and dCOGs
AIS_with_diffs = diffs_computer(AIS_with_stats)

### Calculating summary statistics ###

# Assigning new pointer, releasing old memory
AIS_db = AIS_with_diffs

# Calculating summary statistics
AIS_summary_stats = statistics_aggregator(AIS_db)

# Finally writing output for further analysis
AIS_summary_stats.to_csv("AIS_summary_stats.csv",index=False,
                     )

