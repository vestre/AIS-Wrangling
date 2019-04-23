"""
File for executing whole AIS workflow
"""

# Importing relevant libraries
import numpy as np
import pandas as pd 
import datetime

from NearCollisionStatistics *

### Reading from file and sorting ###

# Setting internal options for speedup
pd.set_option('compute.use_bottleneck', True)
pd.set_option('compute.use_numexpr', True)

PATH = "~/DATA/aisdata.csv"

# Setting columns to write
incols = [1,2,3,4,5,6,7,8,17]
incols2 = ["date_time_utc","mmsi","lat","lon","sog","cog","true_heading","length","nav_status"]

# Setting datecolumn if we want to convert times directly
datecol = 2

# Reading file
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
                                "nav_status": np.int32}
                        )

# Renaming columns
AIS_df.rename(columns={"date_time_utc": "Time",
                          "mmsi":        "ID",
                          "lon":        "LON",
                          "lat":        "LAT",
                          "nav_status": "Status",
                          "true_heading": "Heading",
                          "sog":        "SOG",
                          "cog":        "COG",
                          "length":     "Length"
                         },inplace=True)


# Converting time column from string to datetime
AIS_df["Time"] = pd.to_datetime(AIS_df.Time,format="%Y-%m-%dT%H:%M:%S.%f",box=False)

# Storing datetimes in the dataframe
AIS_df["Time_DateTime"] = AIS_df["Time"]

# Converting datetimes to ints for later conversion to intervals
AIS_df["Time"] = AIS_df.Time.values.astype(np.uint64) / 10**9

# Removing platforms
AIS_df = AIS_df[(AIS_df.ID > 99999999)&(AIS_df.ID <= 999999999)]

# Sorting by time
AIS_df.sort_values(by = "Time",inplace=True)

# Writing file to HDF
AIS_df.to_hdf("AIS_NOR.hdf",
              "AIS",
              mode = "w",
              format="table")

### Finding unique vessels ###

PATH_hdf1 = "/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/computation/AIS_NOR.hdf"

AIS_sorted = pd.read_hdf(PATH_hdf1,key="AIS")

# Creating intervals of 20 seconds
AIS_sorted["Interval"] = (AIS_sorted.Time - np.min(AIS_sorted.Time)) // 20

# Option: Interval of 60 minutes
AIS_sorted["IntervalMins"] = (AIS_sorted.Time - np.min(AIS_sorted.Time)) // 60

# Making Interval iterable (to integer values)
AIS_sorted.Interval = AIS_sorted.Interval.values.astype(np.uint32)
AIS_sorted.IntervalMins = AIS_sorted.IntervalMins.values.astype(np.uint32)

# Dropping duplicates if they have both the same ID and the same intervalMin, 
# minus the first instance - going from 18 to 8 M instances
AIS_shorter = AIS_sorted.drop_duplicates(["ID","Interval"],keep="first")

# Sorting in order to speed up computations (making LAT and LON ordinal)
AIS_shorter = AIS_shorter.sort_values(by=["Time","LAT","LON"], axis=0, inplace=False)

# Saving the new file
AIS_shorter.to_hdf("AIS_NOR_unique_20sec.hdf",
              "AIS",
              mode = "w",
              format="table")

### Computing times to CPA and CPA ###

PATH_hdf2 = "/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/computation/AIS_NOR_unique_20sec.hdf"

AIS_shorter = pd.read_hdf(PATH_hdf2,key="AIS")

AIS_outmatrix = time_to_CPA_calculator(AIS_shorter,extrainterval=False)

# For internal storage
AIS_outmatrix.to_hdf("Filtered_AIS_tCPA_20sec_improved.hdf",
              "AIS_candidates",
              mode = "w",
              format="table")

### Filtering the data by tCPA and CPA ###

# Re-reading files
PATH_csv = "/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/computation/Filtered_AIS_tCPA_20sec_improved.hdf"
AIS_filtered = pd.read_hdf(PATH_csv)

# Filtering by SOG and Time to CPA and Navigation status
SOG_tresh_min = 5
SOG_tresh_max = 50
filter3 = (((AIS_outmatrix["SOG_Vessel_1"]>SOG_tresh_min)&(AIS_outmatrix["SOG_Vessel_2"]>SOG_tresh_min))
            &((AIS_outmatrix["SOG_Vessel_1"]<SOG_tresh_max)&(AIS_outmatrix["SOG_Vessel_2"]<SOG_tresh_max))
            &((AIS_outmatrix["Status_Vessel_1"] != 11)&(AIS_outmatrix["Status_Vessel_1"] != 12))
            &((AIS_outmatrix["Status_Vessel_2"] != 11)&(AIS_outmatrix["Status_Vessel_2"] != 12))
            &(AIS_outmatrix["Min_time_to_CPA"] < 1000)&(AIS_outmatrix["CPA"] < 1000))


AIS_filtered_more = AIS_outmatrix[filter3]

# Sorting and finding unique situations
AIS_Filtered_more = AIS_filtered_more.sort_values(["Min_time_to_CPA"])
AIS_Unique_Filtered_more = AIS_Filtered_more.drop_duplicates(["ID_Vessel_1","ID_Vessel_2"])
AIS_Unique_Filtered_more = AIS_Unique_Filtered_more.sort_values(["Min_time_to_CPA"])

# For internal storage
AIS_Unique_Filtered_more.to_hdf("Filtered_AIS_20s_tight.hdf",
              "AIS_candidates",
              mode = "w",
              format="table")


### Building a databank ###

# Re-reading unique filtered situation data
PATH_hdf = "/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/computation/Final_out/Filtered_AIS_20s_tight.hdf"
AIS_filtered = pd.read_hdf(PATH_hdf, key = "AIS_candidates")

# Re-reading full AIS original file
PATH = "~/code/DATA_DNV/ais_small.csv"

incols = ["date_time_utc","mmsi","lat","lon","sog","cog","true_heading","length","nav_status","RISK_Norwegian_Main_Vessel_Category_ID"]
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
AIS_df["Time"] = pd.to_datetime(AIS_df.Time,format="%Y-%m-%dT%H:%M:%S.%f",box=False)

# Storing datetimes in the dataframe
AIS_df["Time_datetime"] = AIS_df["Time"]

# Converting datetimes to ints for later conversion to intervals
AIS_df["Time"] = AIS_df.Time.values.astype(np.uint64) / 10**9

# Removing platforms
AIS_df = AIS_df[(AIS_df.ID > 99999999)&(AIS_df.ID <= 999999999)]

# Sorting by time
AIS_df.sort_values(by = "Time",inplace=True)

# Collecting observations
DATABANK = observation_collector(AIS_df, AIS_filtered)

# Writing to file
DATABANK.to_hdf("AIS_Databank_tight.hdf",
              "AIS",
              mode = "w",
              format="table")

### Synchronizing times ###

# Re-reading files
PATH_hdf = "/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/computation/RecordsAnalysis/AIS_Databank_tight.hdf"
AIS_databank = pd.read_hdf(PATH_hdf, key = "AIS")

# Synchronizing times
databank_mergetime = observation_synchronizer(AIS_databank)

# Writing checkpoint
databank_mergetime.to_hdf("AIS_databank_mergedtimes_comb.hdf",
              "AIS",
              mode = "w",
              format="table")

### Calculating Statistics (CPA, tCPA, CPA_dist) for all observations ###

# Re-reading databank
PATH_hdf = "/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/computation/RecordsAnalysis/AIS_databank_mergedtimes_comb.hdf"
AIS_db = pd.read_hdf(PATH_hdf, key = "AIS")

# Calculating stats
AIS_with_stats = Stat_calculator(AIS_db)

# Calculating dSOG and dCOGs
AIS_with_diffs = diffs_computer(AIS_with_stats)

# Saving checkpoint
AIS_with_diffs.to_hdf("AIS_with_stats.hdf",
              "AIS",
              mode = "w",
              format="table")

### Calculating summary statistics ###

# Re-reading files
PATH_hdf = "/Users/arnsteinvestre/Desktop/Studier/19_vår/STK-MAT2011/3_Files/FINAL_NOR/computation/RecordsAnalysis/AIS_with_stats.hdf"
AIS_db = pd.read_hdf(PATH_hdf, key = "AIS")

# Calculating summart statistics
AIS_summary_stats = statistics_aggregator(AIS_db)

# Writing final output for further analysis
AIS_summary_stats.to_csv("AIS_summary_stats.csv",index=False,
                     )

