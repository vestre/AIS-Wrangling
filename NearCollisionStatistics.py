"""
Module for making near-collision situation summary statistics.
"""

import numpy as np 
import pandas as pd
import numba

pd.set_option('compute.use_bottleneck', True)
pd.set_option('compute.use_numexpr', True)

def compute_nearcollision_statistics(
                PATH,
                incols = ["date_time_utc","mmsi","lat","lon","sog","cog",
                                "true_heading","length","nav_status"], 
                namedict = {"DateTime": "Time_Datetime", 
                            "ID": "ID",
                            "Time": "Time",
                            "COG": "COG",
                            "SOG": "SOG",
                            "Heading": "Heading",
                            "LON": "LON",
                            "LAT": "LAT",
                            "Status": "Status",
                            "Length": "Length",
                            "Interval": "Interval",
                            "Distances": "Distances"},
                colnamer={"date_time_utc": "Time",
                              "mmsi":        "ID",
                              "lon":        "LON",
                              "lat":        "LAT",
                              "nav_status": "Status",
                              "true_heading": "Heading",
                              "sog":        "SOG",
                              "cog":        "COG",
                              "length":     "Length"
                             },
                dtypesforreading={"mmsi": np.uint64,             
                                "lat": np.float32, 
                                "lon": np.float32,
                                "sog": np.float32,
                                "cog": np.float32,
                                "true_heading": np.float32,
                                "length": np.float32,
                                "nav_status": np.int32}
                intervalcol=False,
                extrainterval=False):
    """
    Takes path to table, name of relevant columns in dict and to read as well
    as further options. Outputs summary statistics.

    NB: Still beta.
    """
    print("Reading csv file from PATH")
    AIS_df = pd.io.parsers.read_csv(PATH
                        ,engine="c"
                        ,sep=";"
                        ,usecols=incols        
                        ,compression=None
                        ,dtype=dtypesforreading,
                        )

    print("csv file read")

    # Renaming columns
    AIS_df.rename(columns=colnamer,inplace=True)

    # Converting time column from string to datetime
    AIS_df[namedict["Time"]] = pd.to_datetime(AIS_df.Time,format="%Y-%m-%dT%H:%M:%S.%f",box=False)

    # Storing datetimes in the dataframe
    AIS_df[namedict["DateTime"]] = AIS_df["Time"]

    # Converting datetimes to ints for later conversion to intervals
    AIS_df[namedict["Time"]] = AIS_df.Time.values.astype(np.uint64) / 10**9

    # Removing platforms
    AIS_df = AIS_df[(AIS_df.ID > 99999999)&(AIS_df.ID <= 999999999)]

    # Sorting by time
    AIS_df.sort_values(by = namedict["Time"],inplace=True)

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
    print("Computing time to CPA for all records")
    # Calculating times to CPA and CPA distances for all values in all
    # intervals
    # NB: This is the most computationally heavy operation - see function code
    AIS_outmatrix = time_to_CPA_calculator(AIS_shorter, extrainterval=False,
                                            namedict = namedict)
    print("All CPAs computed")
    print("Filtering outtable")
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
    print("Table filtered. Creating databank.")

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

    print("Synchronizing databank")
    # Synchronizing times
    databank_mergetime = observation_synchronizer(AIS_databank)

    ### Calculating Statistics (CPA, tCPA, CPA_dist) for all observations ###

    # Creating pointer
    AIS_db = databank_mergetime

    print("Calculating statistics")

    # Calculating stats
    AIS_with_stats = Stat_calculator(AIS_db)

    # Calculating dSOG and dCOGs
    AIS_with_diffs = diffs_computer(AIS_with_stats)

    ### Calculating summary statistics ###

    # Assigning new pointer, releasing old memory
    AIS_db = AIS_with_diffs

    # Calculating summary statistics
    AIS_summary_stats = statistics_aggregator(AIS_db)

    print("Writing statistics to file.")

    # Finally writing output for further analysis
    AIS_summary_stats.to_csv("AIS_summary_stats.csv",index=False,
                         )

def time_to_CPA_calculator(whole_table, 
                namedict = {"DateTime": "Time_Datetime", 
                            "ID": "ID",
                            "Time": "Time",
                            "COG": "COG",
                            "SOG": "SOG",
                            "Heading": "Heading",
                            "LON": "LON",
                            "LAT": "LAT",
                            "Status": "Status",
                            "Length": "Length",
                            "Interval": "Interval",
                            "Distances": "Distances"},
                intervalcol=False,
                extrainterval=False):
    """
    Wrapper function for Time to CPA calculation. Takes pandas
    dataframe (table) and a list of column names, in order to
    generalize the extraction of table series as arrays to
    the workhorse calculation function.

    The interval field is made by integer division of Time (in seconds)
    by selected interval cutoff (e.g. 20 seconds) pre application of
    the function.

    Function works well with dataframe where unique instances / records
    are found for every ID in every interval (by pd.df.drop_duplicates), 
    outside of the function.

    Input: Pandas dataframe with records, dictionary translating
    column names to pre-specified column names used below in the 
    function

    Output: Pandas dataframe with Ship IDs, Timestamps, ID to ship
    with which the ship has its closest time to CPA, time to CPA at 
    smallest value and distance at this closest time value.
    """

    # Optional argument to set the interval column name directly to
    # the function. Useful for selecting other interval sizes if these
    # are defined.
    if intervalcol == False:
        pass
    else:
        namedict["Interval"] = intervalcol

    # Setting minimum tolerance for close-to-zero bool tests
    # and meters to degree conversion constant
    TOL = 1E-7
    METERS_DEGREES_CONVERSION = 1.1132e+5

    # Make sure to sort table by interval, LON and LAT. Sorting by LON as
    # second enables us to break loops further down
    whole_table = whole_table.sort_values(
                                by=[namedict["Interval"],
                                    namedict["LON"],
                                    namedict["LAT"]], 
                                axis=0, 
                                inplace=False)

    # Making set of indices for later use in extracting the correct
    # information for ship 2
    Indices = np.arange(len(whole_table.index))
    whole_table.set_index(Indices,inplace=True)

    # Calls workhorse function
    outmatrix = _CPA_calculator_workhorse(
        np.asarray(whole_table[namedict["ID"]]),
        np.asarray(whole_table[namedict["LAT"]]),
        np.asarray(whole_table[namedict["LON"]]),
        np.asarray(whole_table[namedict["Time"]]),
        np.asarray(whole_table[namedict["SOG"]]),
        np.asarray(whole_table[namedict["COG"]])*np.pi/180,
        np.asarray(whole_table[namedict["Interval"]]),
        np.asarray(Indices),
        extrainterval=extrainterval
        )

    # Creates mask in order to filter out all zero or close-to-zero
    # values of time to CPA and negative values of CPA. 
    # Use this to truncate outmatrix to only get interesting values 
    # (zero is used as "pass" value in the workhorse function)
    zeros_mask = (outmatrix[2] > TOL)&(outmatrix[3] >= 0)
    indices_mins = outmatrix[1]

    # Collecting pointers to the series from the main dataframe
    IDs = np.asarray(whole_table[namedict["ID"]])
    Times  = np.asarray(whole_table[namedict["DateTime"]])
    COGs = np.asarray(whole_table[namedict["COG"]])
    SOGs = np.asarray(whole_table[namedict["SOG"]])
    Status = np.asarray(whole_table[namedict["Status"]])

    # Collecting data for ship 1. These are collected directly by using
    # that we only want the zero mask-data (tCPA > TOL and CPA>=0).
    ID_V1 = outmatrix[0][zeros_mask]
    Time_V1 = Times[zeros_mask]
    COG_V1 = COGs[zeros_mask]
    SOG_V1 = SOGs[zeros_mask]
    Status_V1 = Status[zeros_mask]

    # Creating arrays for the closest ship 2 instances by using 
    # the indices collected with the workhorse function
    ID_V2_tot = IDs[indices_mins]
    Time_V2_tot = Times[indices_mins]
    COG_V2_tot = COGs[indices_mins]
    SOG_V2_tot = SOGs[indices_mins]
    Status_V2_tot = Status[indices_mins]

    # Collecting the final subset of data for ship 2 by removing all
    # instances where the tCPA and CPA distance indicate no
    # danger of collision
    ID_V2 = ID_V2_tot[zeros_mask]
    Time_V2 = Time_V2_tot[zeros_mask]
    COG_V2 = COG_V2_tot[zeros_mask]
    SOG_V2 = SOG_V2_tot[zeros_mask]
    Status_V2 = Status_V2_tot[zeros_mask]

    # Finally collecting tCPA and CPA, making sure to convert CPA_dist from
    # degrees to meters and filtering only interesting instances as above
    Time_to_CPA = outmatrix[2][zeros_mask]
    CPA_dist = outmatrix[3][zeros_mask]*METERS_DEGREES_CONVERSION

    # Creating output object as dataframe
    outdataframe = pd.DataFrame(data = 
                                {"ID_Vessel_1":     ID_V1,
                                "Time_Vessel_1":    Time_V1,
                                "SOG_Vessel_1":     SOG_V1,
                                "COG_Vessel_1":     COG_V1,
                                "Status_Vessel_1":  Status_V1,
                                "ID_Vessel_2":      ID_V2,
                                "Time_Vessel_2":    Time_V2,
                                "SOG_Vessel_2":     SOG_V2,
                                "COG_Vessel_2":     COG_V2,
                                "Status_Vessel_2":  Status_V2,
                                "Min_time_to_CPA":  Time_to_CPA,
                                "CPA":              CPA_dist}, 
                        columns = ["ID_Vessel_1",
                                "Time_Vessel_1",
                                "SOG_Vessel_1" ,
                                "COG_Vessel_1" ,
                                "Status_Vessel_1" ,
                                "ID_Vessel_2" ,
                                "Time_Vessel_2",
                                "SOG_Vessel_2" ,
                                "COG_Vessel_2" ,
                                "Status_Vessel_2" ,
                                "Min_time_to_CPA",
                                "CPA"])

    return outdataframe

@numba.jit(parallel=True, nopython=True)
def _CPA_calculator_workhorse(IDs, LATs, LONs, DateTimes,
                            SOGs, COGs, Intervals, Indices,
                            extrainterval=False):
    """
    Workhorse function for time to closest point of approach-computation.
    Takes lists of corresponding records of ID, timestamps, latitudes, 
    longitudes, speed over ground, course over ground for ships, as well as an 
    index separating the data into time intervals of a certain length, as an 
    array with the distance from each ship to the closest ship. Returns a matrix
    of ship IDs, IDs of the other ship with which the ship has the smallest
    time to closest point of approach, the time to closest point of approach and
    the closest point of approach. Only computes for relevant collision
    candidates

    Input: Arrays of IDs, Latitudes, Longitudes, speed over ground, course
    over ground, interval index, dataframe index and optional argument for 
    double interval iteration

    Output: Matrix with IDs, index of ship with closest tCPA, time to CPA, 
    CPA-distance for all records in original frame
    """

    # Setting minimum distance to 2000 M, and converting
    # minimum distance to meters. NB: Will be latitude-
    # corrected at later stage
    METERS_DEGREES_CONVERSION = 1.1132e+5
    DIST_MIN = 2000/METERS_DEGREES_CONVERSION

    # Setting minimum tolerance for close-to-zero bool tests
    TOL = 1E-7

    # Creates index numbers for interval parsing and save length of
    # table for array creation
    interval_numbers = np.unique(Intervals)
    first_interval = interval_numbers[0]
    table_length = len(IDs)

    # Create arrays for collection of results
    CPA_dist = np.zeros(table_length, dtype=np.float32)
    CPA_time = np.zeros(table_length, dtype=np.float32)
    Indices_mins = np.zeros(table_length, dtype=np.uint32)

    # Iterating over intervals for interval-wise calculations
    # (written to be implemented with parallellism)
    for interval in numba.prange(first_interval,
                        first_interval + len(interval_numbers)):

        # Calculating interval size
        interval_start = np.uint32(np.sum(Intervals < interval))
        interval_size = np.uint32(np.sum(Intervals == interval))

        if extrainterval == True:
            interval_next_size = np.uint32(np.sum(Intervals == interval + 1))
        else:
            interval_next_size = 0

        long_interval_size = interval_size + interval_next_size

        for record_index1 in range(interval_start,interval_start 
                                                    + interval_size):

            # If ship is not already checked, we establish speed, angle,
            # LON and LAT for ship 1
            ship1_lon = LONs[record_index1]
            ship1_lat = LATs[record_index1]
            ship1_pos = np.array([ship1_lon, ship1_lat],
                                    dtype=np.float32)
            speed1 = SOGs[record_index1]
            angle1 = COGs[record_index1]

            # Creating arrays for storing CPA_distance and time to
            # CPA for 1he iteration over ship 1
            CPA_ship_1 = np.zeros(long_interval_size, dtype = np.float32)
            tCPA_ship_1 = np.zeros(long_interval_size, dtype = np.float32)

            # Iterating from ship 1 to end of interval to calculate
            # time to CPA to all other (relevant) ships
            for record_index2 in range(record_index1+1, interval_start 
                                                        + long_interval_size):

                # Using the fact that the table is sorted by Longitudes to
                # break the loop if we have passed by the relevant LON area.
                # Passing smaller lons
                if (LONs[record_index2] > ship1_lon + DIST_MIN):
                    break
                elif (LONs[record_index2] < ship1_lon - DIST_MIN):
                    pass
                else:
                    # Creating interval index for ship 2
                    # Establishing LON and LAT for ship 2
                    ship2_lon = LONs[record_index2]
                    ship2_lat = LATs[record_index2]

                    # Establishing LAT and LON mean value for correct
                    # latitude correction value (computed in radians)
                    LAT_mean = np.pi*(ship2_lat + ship1_lat)/(2*180)

                    # Correcting minimum distance for 
                    DIST_MIN_CORR = DIST_MIN/np.cos(LAT_mean)

                    # We check wether the two ships are close in LAT
                    if ship2_lat > ship1_lat + DIST_MIN_CORR:
                        pass

                    elif ship2_lat < ship1_lat - DIST_MIN_CORR:
                        pass

                    else:
                        # Making initial position vectors, in degrees
                        ship2_pos = np.array([ship2_lon, ship2_lat],
                                                dtype=np.float32)

                        # Making initial velocity vectors, measuring in degrees
                        # per seconds (converting from "speed" in knots)
                        velocity1 = (speed1
                                        *np.array([np.sin(angle1),
                                          np.cos(angle1)],
                                                    dtype=np.float32))

                        velocity2 = (SOGs[record_index2]
                                        *np.array([np.sin(COGs[record_index2]),
                                          np.cos(COGs[record_index2])],
                                                    dtype=np.float32))

                        # Calculating initial distance and crossing velocities
                        # (evolution in distance between ships in plane)
                        distance_at_zero = ship1_pos - ship2_pos
                        cross_velocities = velocity1 - velocity2

                        # Correcting distances and velocities in LON direction 
                        # for the mean latitude skewness 
                        distance_at_zero[0] = (distance_at_zero[0]
                                                *np.cos(LAT_mean))
                        cross_velocities[0] = (cross_velocities[0]
                                                *np.cos(LAT_mean))

                        # Calculating the squared norm of the CV for future use
                        cross_velocity_squared_norm = (
                                 (cross_velocities[0])*(cross_velocities[0])
                               + (cross_velocities[1])*(cross_velocities[1]))

                        # Setting time_to_CPA as 0 if ships have identical 
                        # velocities or negative evolution in distance in 
                        # both directions. Calculating Time to CPA if not.
                        if (cross_velocity_squared_norm < TOL):
                            time_to_CPA = 0

                        else:
                            time_to_CPA = -(
                                (distance_at_zero[0]*cross_velocities[0]
                                + distance_at_zero[1]*cross_velocities[1])
                                / cross_velocity_squared_norm)

                            # Converting to seconds by formula from
                            # Azzeddine
                            time_to_CPA_sec = (time_to_CPA 
                                                *METERS_DEGREES_CONVERSION 
                                                / 0.514444)

                        # Setting CPA negative if ships are heading in opposite 
                        # directions (time_to_CPA is negative), zero or above 
                        # 20 minutes. Calculating CPA if not.
                        if (time_to_CPA_sec < TOL or time_to_CPA_sec > (40*60)):
                            distance_to_CPA = -1.0

                        else:
                            distvec = (distance_at_zero 
                                    + time_to_CPA*cross_velocities)

                            distance_to_CPA = np.sqrt(distvec[0]*distvec[0]
                                                    + distvec[1]*distvec[1])

                        # Saving relevant values in ship 1 collection array
                        tCPA_ship_1[record_index2 - interval_start] = (
                                                    np.float32(time_to_CPA_sec))
                        CPA_ship_1[record_index2 - interval_start] = ( 
                                                    np.float32(distance_to_CPA))

            # Creating filter for whether ship 1 has any relevant time_to_CPA
            # situations
            tCPA_relevant_filter = tCPA_ship_1 > 0

            # Using filter to only save if there are relevant situtations
            if np.any(tCPA_relevant_filter):

                # Creating candidates with the filters so that the 
                # we are not taking the minimum of zero or negative values
                # and so that min_tCPA_ind will work as index in CPA_dist
                min_tCPA_candidates = tCPA_ship_1[tCPA_relevant_filter]
                min_CPAdist_candidates = CPA_ship_1[tCPA_relevant_filter]

                # Creating indices list in order to be able to get correct
                # index number by min_tCPA_ind
                Indices_candidates = Indices[Intervals == interval]
                Indices_candidates = Indices_candidates[tCPA_relevant_filter]

                # Finding the index for the minimum tCPA
                min_tCPA_ind = np.argmin(min_tCPA_candidates)

                # Collecting relevant data for the smallest-time-to-CPA 
                # incidence in the global function collection arrays
                CPA_time[record_index1] = np.float32(
                                        min_tCPA_candidates[min_tCPA_ind])
                
                CPA_dist[record_index1] = np.float32(
                                        min_CPAdist_candidates[min_tCPA_ind])
                
                Indices_mins[record_index1] = np.uint32(
                                        Indices_candidates[min_tCPA_ind])


    return (IDs,
            Indices_mins,
            CPA_time, 
            CPA_dist)


def observation_collector(bigtable, selectortable,
                            idcolsel1 = "ID_Vessel_1", 
                            idcolsel2 = "ID_Vessel_2",  
                            timecolsel1 = "Time_Vessel_1",  
                            timecolsel2 = "Time_Vessel_2",
                            bigtableIDname = "ID",
                            bigtableTimename = "Time",
                            bigtableCatname = "Category",
                            timewindow = 3600): 

    """
    Retrieves near-collision situations from large in-table based on list of
    near-collision situation in selection table, previously determined.

    Input: Large AIS data table, selection table, column name 
    specifiers (optional)

    Output: New table of records for near-collision situations.
    """


    # Setting one side of the time window with pre-defined value (set at 1 hour
    # as standard)
    timewindow = timewindow

    # Excluding IDs not in selector and wrong category
    filter_bigtab = ((np.isin(bigtable[bigtableIDname],selectortable[idcolsel1])
                |np.isin(bigtable[bigtableIDname],selectortable[idcolsel2]))
                &(bigtable["Category"] != 12)&(bigtable["Category"].notna()))

    # Removing uninteresting records
    bigtable = bigtable[filter_bigtab]

    # Making filter for the relevant IDs for selectortable 
    filter_selec = (np.isin(selectortable[idcolsel1],bigtable[bigtableIDname])
                    &np.isin(selectortable[idcolsel2],bigtable[bigtableIDname]))

    # Removing uninteresting records 
    selectortable = selectortable[filter_selec]

    # Sorting by time
    bigtable = bigtable.sort_values(["Time"])

    # Resetting indexing to be sure of correct picking of 
    # content at end
    selectortable = selectortable.reset_index(drop=True)
    bigtable = bigtable.reset_index(drop=True) 

    bigtable_index = bigtable.index

    # Collecting boolean filter vector, and labels for later filtering
    # NB: Bigtable[Time] is already converted to int64. That is not the case
    # in the selector table
    fullselector, sitlist = _index_finder_workhorse(
                        np.asarray(selectortable[idcolsel1]),
                        np.asarray(selectortable[idcolsel2]),
                        np.asarray(selectortable[timecolsel1]
                            .astype(np.int64)/10**9).astype(np.int64),
                        np.asarray(selectortable[timecolsel2]
                            .astype(np.int64)/10**9).astype(np.int64),
                        np.asarray(bigtable[bigtableIDname]),
                        np.asarray(bigtable[bigtableTimename]).astype(np.int64),
                        bigtablen = len(bigtable.index),
                        seltablen = len(selectortable.index),
                        timewindow = timewindow,
                        bigtable_index = np.asarray(bigtable_index))

    # Filtering table
    outtable = bigtable.iloc[fullselector]

    # Making shallow copy for making new columns
    outtable = outtable.copy(deep=False)

    # Storing labels for situtation number and ship in encounter
    outtable["Situations"] = sitlist

    # Sorting table for pretty output
    outtable = outtable.sort_values(["Situations", "ID", "Time"])

    return outtable

@numba.jit(nopython=False)
def _index_finder_workhorse(IDs1, IDs2,
                            TimeSel1, TimeSel2,
                            bigtableID,
                            bigtableTime,
                            bigtablen,
                            seltablen,
                            timewindow,
                            bigtable_index):
    """
    Workhorse function for determining the correct indices for the situation
    and outputing these for collection from the full table.

    NB: Not fully optimized, but finishes in OK time.
    """

    # Creating selectors and collectors
    fullselector = np.zeros(0, dtype = np.uint32)
    sitlist = np.zeros(0, dtype = np.uint32)
    shiplist = np.zeros(0, dtype = np.uint32)

    # Iterating over selector table and withdrawing close instances from
    # Full table
    for i in range(seltablen):

        # Selecting relevant ID and time points to store indices
        selector1 = ((bigtableID == IDs1[i])
            &(np.abs(bigtableTime - TimeSel1[i]) < timewindow))

        selector2 = ((bigtableID == IDs2[i])
            &(np.abs(bigtableTime - TimeSel2[i]) < timewindow))

        # Aggregating boolean indices
        fullselector = np.concatenate([fullselector , 
                            bigtable_index[selector1|selector2]], axis = 0)

        # Storing situation number
        sitlist = np.concatenate([sitlist, np.uint32(i)*np.ones(
                            np.sum(selector1|selector2),np.uint32)], axis = 0)

    return fullselector, sitlist



def observation_synchronizer(whole_table, 
                            timecolname = "Time", 
                            sitcolname = "Situations", 
                            idcolname = "ID",
                            statcolname = "Status",
                            SOGcolname = "SOG",
                            COGcolname = "COG",
                            catcolname = "Category",
                            datetimcol = "Time_datetime",
                            LATcolname = "LAT",
                            LONcolname = "LON",
                            lencolname = "Length",
                            Headcolname = "Heading"
                            ):
    """
    Function for synchronizing observations in near-collision avoidance 
    situations where records are already collected. Also combines situations
    so that the output give each point in time in a situation as one row with
    information for both vessels.

    Input: Table of records for each vessel in each situation, labelled by ID.
    Output: Table of record-pairs for each point in time, with synced time
    scale.
    """

    # Create collector tables for situation tables
    tablecollector1 = []
    tablecollector2 = []

    # Looping over all unique siutations
    for i in np.unique(whole_table[sitcolname]):

        # Retrieving all records pertaining to each one situation
        situation_table = whole_table[whole_table[sitcolname] == i]

        # Determining the IDs of the vessels in the situation
        IDs = np.unique(situation_table[idcolname])

        # Asserting the tables of records for each of the vessels in the 
        # situation
        table1 = situation_table[situation_table[idcolname] == IDs[0]]
        table2 = situation_table[situation_table[idcolname] == IDs[1]]

        # Resetting indices for each of the two tables, to ensure no indexing
        # trouble
        table1 = table1.reset_index(drop=True)
        table2 = table2.reset_index(drop=True)

        # Determining the shorter and longer time series of records
        if len(table1.index) <= len(table2.index):
            shortest_table_index = table1.index
            longest_table_index = table2.index
        else:
            shortest_table_index = table2.index
            longest_table_index = table1.index

        # Synchronizing times with workhorse function
        short_table_ind, long_table_ind, mergetimes = _obs_sync_workhorse(
                            np.asarray(table1[timecolname]), 
                            np.asarray(table2[timecolname]), 
                            time1_index = np.asarray(shortest_table_index),
                            time1_len = len(shortest_table_index),
                            time2_len = len(longest_table_index))

        # Collecting by the correct long/short table
        if len(table1.index) <= len(table2.index):
            table1 = table1.iloc[short_table_ind]
            table2 = table2.iloc[long_table_ind]
        else:
            table1 = table1.iloc[long_table_ind]
            table2 = table2.iloc[short_table_ind]

        # Saving the new merged time scale
        table1["Merge_times"] = mergetimes

        # Adding the tables to the two collectors for future concatenation
        tablecollector1.append(table1)
        tablecollector2.append(table2)

    # Concatenating tables for each of the two sets of situations
    out1 = pd.concat(tablecollector1, ignore_index=True, copy=False, sort=False,axis=0)
    out2 = pd.concat(tablecollector2, ignore_index=True, copy=False, sort=False,axis=0)

    # Renaming columns for both concatenated tables
    out1.rename({timecolname: "Time_1",
                    idcolname: "ID_1",
                    statcolname: "Status_1",
                    COGcolname: "COG_1",
                    SOGcolname: "SOG_1",
                    catcolname: "Category_1",
                    LATcolname: "LAT_1",
                    LONcolname: "LON_1",
                    Headcolname: "Heading_1",
                    lencolname: "Length_1",
                    sitcolname: "Situations_1",
                    datetimcol: "Time_datetime_1"
                    },axis="columns",
                    inplace=True)

    out2.rename({timecolname: "Time_2",
                    idcolname: "ID_2",
                    statcolname: "Status_2",
                    COGcolname: "COG_2",
                    SOGcolname: "SOG_2",
                    catcolname: "Category_2",
                    LATcolname: "LAT_2",
                    LONcolname: "LON_2",
                    Headcolname: "Heading_2",
                    lencolname: "Length_2",
                    sitcolname: "Situations_2",
                    datetimcol: "Time_datetime_2"
                    },axis="columns",
                    inplace=True)

    # Concatenating the two "collected" tables, but this time by columns
    outtable = pd.concat([out1, out2], 
                            axis=1, 
                            ignore_index=False, 
                            copy=False, 
                            sort=False)

    # Asserting the datetime format of merge times (remembering to multiply
    # with 10**9 in order to take into account the nanosecond format of
    # pandas datatime)
    outtable["MT_datetime"] = pd.to_datetime(outtable["Merge_times"]*10**9)

    # Setting the order in the new output merged table
    newcolorder = ["ID_1","ID_2","MT_datetime","LAT_1","LAT_2",
                    "LON_1","LON_2","COG_1","Heading_1","COG_2",
                    "Time_1","Time_2","Time_datetime_1",
                    "Time_datetime_2","Merge_times","Heading_2",
                    "SOG_1","SOG_2","Status_1","Status_2",
                    "Category_1","Category_2","Length_1","Length_2",
                    "Situations_1","Situations_2"]

    # Reindexing by the new column order
    outtable = outtable.reindex(newcolorder, axis = "columns",copy=False)

    # Changing types of the integer values (this has changed somewhere
    # in the process)
    outtable["ID_1"] = outtable["ID_1"].astype(np.uint32)
    outtable["ID_2"] = outtable["ID_2"].astype(np.uint32)
    outtable["Status_1"] = outtable["Status_1"].astype(np.uint32)
    outtable["Status_2"] = outtable["Status_2"].astype(np.uint32)
    outtable["Situations_1"] = outtable["Situations_1"].astype(np.uint32)   
    outtable["Situations_2"] = outtable["Situations_2"].astype(np.uint32)   

    return outtable

@numba.jit(nopython=True, parallel=False)
def _obs_sync_workhorse(intimes1, intimes2, time1_index, time1_len, time2_len):
    """
    Workhorse function that takes time values and produces a synchronized set of
    indices for the two in-tables, as well as a mean time array.
    """
    # Sets table1 as the smallest table (for the purpose of this sub-function)
    if time1_len <= time2_len:
        time1 = np.asarray(intimes1,dtype = np.int64)
        time2 = np.asarray(intimes2,dtype = np.int64)
    else:
        time1 = np.asarray(intimes2,dtype = np.int64)
        time2 = np.asarray(intimes1,dtype = np.int64)

    # Asserting collection array in order to filter out those time values in 
    # the shorter table where there is no corresponding value in the longer
    # table within plus-minus 5 seconds
    notime_filter = np.ones(time1_len, dtype = np.bool_)

    # Asserting collection array for indices for table 2 to fund the closest
    # times. These will continue to be zero for all times which are not
    # found to be close
    closetimes_ind = np.zeros(time2_len,dtype = np.int32) - 1

    # Looping through the smallest time array in order to match each point
    # with the closest point in the longer array
    for i in numba.prange(time1_len):
        # Finding the index for the closest time point
        ctime_ind = np.abs(time2 - time1[i]).argmin()

        # Filtering out if the closest time is not close enough (5 secs 
        # plus-minus)
        if time1[i] < time2[ctime_ind] - 5:
            notime_filter[i] = False
        elif time1[i] > time2[ctime_ind] + 5:
            notime_filter[i] = False

        # If the situation is not filtered out, the index is instead saved
        else:
            closetimes_ind[i] = ctime_ind

    # Deleting all times with no close counterpart
    table1_index = time1_index[notime_filter]

    # Retrieving all close time point indices from the long table
    table2_index = closetimes_ind[closetimes_ind >= 0]

    # Collecting the times in order to calculate the merged time
    time1_cand = time1.take(table1_index)
    closetimes = time2.take(table2_index)

    # Calculating the merged time
    mergetimes = (time1_cand + closetimes)/2

    return table1_index, table2_index, mergetimes



def diffs_computer(whole_table, 
                    timemergedcolname = "Merge_times",
                    time1colname = "Time_1", 
                    time2colname = "Time_2", 
                    sitcolname = "Situations_1", 
                    idcolname1 = "ID_1",
                    idcolname2 = "ID_2",
                    SOGcolname1 = "SOG_1",
                    SOGcolname2 = "SOG_2",
                    COGcolname1 = "COG_1",
                    COGcolname2 = "COG_2",
                    LONcolname1 = "LON_1",
                    LONcolname2 = "LON_2",
                    LATcolname1 = "LAT_1",
                    LATcolname2 = "LAT_2"):
    """
    Function for calculating the time derivatives of SOG and COG for each
    vessel in a pair of vessels for each situation in a table of near-
    collision situations. Outputs a table where dSOG and dCOG are new columns, 
    for each vessel.

    Input: Table and colnames
    Output: Table with dSOG and dCOG
    """


    # Filtering for erroneous values
    filter1 = ((whole_table[COGcolname1] <= 360)
                &(whole_table[SOGcolname1] < 50))

    filter2 = ((whole_table[COGcolname2] <= 360)
                &(whole_table[SOGcolname2] < 50))

    filter3 = ((whole_table[LATcolname1] > 54)
                &(whole_table[LONcolname1] > -15)
                &(whole_table[LONcolname1] < 40))

    filter4 = ((whole_table[LATcolname2] > 54)
                &(whole_table[LONcolname2] > -15)
                &(whole_table[LONcolname2] < 40))

    print("Removing", np.sum(~(filter1&filter2&filter3&filter4)),
                                                "erroneous records")

    whole_table = whole_table[filter1&filter2&filter3&filter4]

    # Sorting values in order for further computational choices to work.
    # Resetting index in order to make sure indexing works as it should
    whole_table = whole_table.sort_values([sitcolname,timemergedcolname])
    whole_table = whole_table.reset_index(drop=True)

    # Calling workhorse function
    dSOG1s, dSOG2s, dCOG1s, dCOG2s = _diffs_workhorse(
        np.asarray(whole_table[idcolname1]),
        np.asarray(whole_table[idcolname2]),
        np.asarray(whole_table[LATcolname1]),
        np.asarray(whole_table[LATcolname2]),
        np.asarray(whole_table[LONcolname1]),
        np.asarray(whole_table[LONcolname2]),
        np.asarray(whole_table[timemergedcolname]),
        np.asarray(whole_table[time1colname]),
        np.asarray(whole_table[time2colname]),
        np.asarray(whole_table[SOGcolname1]),
        np.asarray(whole_table[SOGcolname2]),
        np.asarray(whole_table[COGcolname1]),
        np.asarray(whole_table[COGcolname2]),
        np.asarray(whole_table[sitcolname])
        )

    # Creating new columns in table with retrieved values
    whole_table["dSOG_1"] = dSOG1s
    whole_table["dSOG_2"] = dSOG2s
    whole_table["dCOG_1"] = dCOG1s
    whole_table["dCOG_2"] = dCOG2s

    # Priting updates on number of records removed
    print("Removing", np.sum((np.isinf(dSOG1s)|np.isinf(dSOG2s)
                        |np.isinf(dCOG1s)|np.isinf(dCOG2s))), "inf entitites")
    print("Removing", np.sum((np.isnan(dSOG1s)|np.isnan(dSOG2s)
                        |np.isnan(dCOG1s)|np.isnan(dCOG2s))), "NaN entitites")

    # Removing records where differentiation has been done by dividing on 0 
    # (creating infs) or 0 on 0 (creating NaNs)
    filtering = ((np.isinf(dSOG1s)|np.isinf(dSOG2s)
                    |np.isinf(dCOG1s)|np.isinf(dCOG2s))
                |(np.isnan(dSOG1s)|np.isnan(dSOG2s)
                    |np.isnan(dCOG1s)|np.isnan(dCOG2s)))

    whole_table = whole_table[~filtering]

    return whole_table

def _diffs_workhorse(ID1s, ID2s, LAT1s, LAT2s, LON1s, LON2s, Mergetimes, 
                Times1, Times2, SOG1s, SOG2s, COG1s, COG2s, Situations):
    """
    Workhorse function that takes table columns as arrays and takes time
    derivatives of the SOGs and GOGs.
    """


    # Creating tolerance for small zero comparisons and creating
    # Conversion constant between degrees and meters
    TOL = 1E-7
    METERS_DEGREES_CONVERSION = 1.1132e+5

    # Asserting table length and creating collector arrays
    table_length = ID1s.size
    dSOG1_collector = np.zeros(table_length,dtype = np.float64)
    dSOG2_collector = np.zeros(table_length,dtype = np.float64)
    dCOG1_collector = np.zeros(table_length,dtype = np.float64)
    dCOG2_collector = np.zeros(table_length,dtype = np.float64)

    # Asserting iteration values for situation loop
    situation_nums = np.unique(Situations)
    first_situation = np.uint32(situation_nums[0])
    number_of_sits = np.uint32(np.max(Situations) + 1)

    for situation in range(first_situation,first_situation + number_of_sits):

        # Asserting iteration values for in-situation loop
        situation_filter = Situations == situation

        # Correcting for whether there is more than two records (enough to
        # differentiate) in the situation
        if np.sum(situation_filter) >= 2:
            with np.errstate(divide='ignore', invalid='ignore'):

                # Differentiating SOGs
                dSOG_1 = (np.diff(SOG1s[situation_filter]) 
                                    / np.diff(Times1[situation_filter]))
                dSOG_2 = (np.diff(SOG2s[situation_filter]) 
                                    / np.diff(Times2[situation_filter]))

                # Setting nans to 0 (they result from 0/0 divison)
                dSOG_1[np.isnan(dSOG_1)] = 0
                dSOG_2[np.isnan(dSOG_2)] = 0

                # Saving dSOGS (need to concatenate with 0 for first value)
                dSOG1_collector[situation_filter] = np.hstack(((0),dSOG_1))
                dSOG2_collector[situation_filter] = np.hstack(((0),dSOG_2))

                # Calculating deltaCOGs (if in (180,360) converted to [0, 180]
                delta_COG1 = np.abs(np.diff(COG1s[situation_filter]) % 360)
                delta_COG1[delta_COG1 > 180]= 360 - delta_COG1[delta_COG1 > 180]

                delta_COG2 = np.abs(np.diff(COG2s[situation_filter]) % 360)
                delta_COG2[delta_COG2 > 180]= 360 - delta_COG2[delta_COG2 > 180]

                # Calculating actual cogs
                dCOG_1 = (delta_COG1 / np.diff(Times1[situation_filter]))
                dCOG_2 = (delta_COG2 / np.diff(Times2[situation_filter]))

                # Same proceedure with Nans and saving
                dCOG_1[np.isnan(dCOG_1)] = 0
                dCOG_2[np.isnan(dCOG_2)] = 0

                dCOG1_collector[situation_filter] = np.hstack(((0),dCOG_1))
                dCOG2_collector[situation_filter] = np.hstack(((0),dCOG_2))

    return (dSOG1_collector, 
            dSOG2_collector, 
            dCOG1_collector, 
            dCOG2_collector)


def Stat_calculator(whole_table, 
                    timemergedcolname = "Merge_times", 
                    sitcolname = "Situations_1", 
                    idcolname1 = "ID_1",
                    idcolname2 = "ID_2",
                    SOGcolname1 = "SOG_1",
                    SOGcolname2 = "SOG_2",
                    COGcolname1 = "COG_1",
                    COGcolname2 = "COG_2",
                    LONcolname1 = "LON_1",
                    LONcolname2 = "LON_2",
                    LATcolname1 = "LAT_1",
                    LATcolname2 = "LAT_2"):
    """
    Wrapper function that takes a table of near-collision situations with ships
    in pair-wise set-up, and a set of pre-set column names (with options). 
    Outputs a table where CPA andtime to CPA is calculated between pairs of 
    vessels

    Input: Table
    Output: Table w CPA and tCPA
    """

    # Sorting values in order for further computational choices to work.
    # Resetting index in order to make sure indexing works as it should
    whole_table = whole_table.sort_values([sitcolname,timemergedcolname])
    whole_table = whole_table.reset_index(drop=True)

    # Calling workhorse function
    (times_to_CPA, CPA_dists,
        dists_to_CPA_1, dists_to_CPA_2
        ) = _stat_workhorse(
        np.asarray(whole_table[idcolname1]),
        np.asarray(whole_table[idcolname2]),
        np.asarray(whole_table[LATcolname1]),
        np.asarray(whole_table[LATcolname2]),
        np.asarray(whole_table[LONcolname1]),
        np.asarray(whole_table[LONcolname2]),
        np.asarray(whole_table[timemergedcolname]),
        np.asarray(whole_table[SOGcolname1]),
        np.asarray(whole_table[SOGcolname2]),
        np.asarray(whole_table[COGcolname1]),
        np.asarray(whole_table[COGcolname2]),
        np.asarray(whole_table[sitcolname])
        )

    # Creating new columns in table with retrieved values
    whole_table["Time_to_CPA"] = times_to_CPA
    whole_table["CPA_distance"] = CPA_dists
    whole_table["Distance_to_CPA_1"] = dists_to_CPA_1
    whole_table["Distance_to_CPA_2"] = dists_to_CPA_2

    return whole_table

@numba.jit(nopython=True, parallel=False)
def _stat_workhorse(ID1s, ID2s, LAT1s, LAT2s, LON1s, LON2s, Mergetimes, 
                SOG1s, SOG2s, COG1s, COG2s, Situations):
    """
    Workhorse function that takes table columns as numpy arrays and makes
    calculations in order to ascertain CPA and time to CPA between every
    pair of vessel in every situation.
    """


    # Creating tolerance for small zero comparisons and creating
    # Conversion constant between degrees and meters
    TOL = 1E-7
    METERS_DEGREES_CONVERSION = 1.1132e+5

    # Asserting table length and creating collector arrays
    table_length = ID1s.size
    time_to_CPA_collector = np.zeros(table_length,dtype = np.float64)
    CPA_dist_collector = np.zeros(table_length,dtype = np.float64)
    dist_to_CPA1_collector = np.zeros(table_length,dtype = np.float64)
    dist_to_CPA2_collector = np.zeros(table_length,dtype = np.float64)

    # Asserting iteration values for situation loop
    situation_nums = np.unique(Situations)
    first_situation = np.uint32(situation_nums[0])
    number_of_sits = np.uint32(np.max(Situations) + 1)

    for situation in range(first_situation, first_situation + number_of_sits):

        # Asserting iteration values for in-situation loop
        situation_start = np.uint32(np.sum(Situations < situation))
        situation_number = np.uint32(np.sum(Situations == situation))

        for i in range(situation_start, situation_start + situation_number):

            # Creating shorthands for indexes in order to simplify notation
            index_ships = i

            ship1_lon = LON1s[index_ships]
            ship1_lat = LAT1s[index_ships]
            ship2_lon = LON2s[index_ships]
            ship2_lat = LAT2s[index_ships]

            speed1 = SOG1s[index_ships]
            speed2 = SOG2s[index_ships]
            angle1 = COG1s[index_ships]
            angle2 = COG2s[index_ships]

            # Calculating mean latitude in order to get distances
            # right in computations
            LAT_mean = np.pi*(ship2_lat + ship1_lat)/(2*180)

            # Asserting velocities as numpy arrays (vectors)
            velocity1 = (speed1
                            *np.array([np.sin(angle1*np.pi/180),
                                        np.cos(angle1*np.pi/180)],
                                        dtype=np.float64))
            velocity2 = (speed2
                            *np.array([np.sin(angle2*np.pi/180),
                                        np.cos(angle2*np.pi/180)],
                                        dtype=np.float64))

            # Asserting distance at point in question
            distance_at_zero = np.array([ship1_lon - ship2_lon, 
                                        ship1_lat - ship2_lat],
                                        dtype=np.float64)

            # Creating cross velocity holder
            cross_velocities = velocity1 - velocity2

            # Correcting distance and veolicy for latitude conversion (LONs
            # are now also skewed)
            distance_at_zero[0] = (distance_at_zero[0]
                                    *np.cos(LAT_mean))
            cross_velocities[0] = (cross_velocities[0]
                                    *np.cos(LAT_mean))

            # Calculating norm for checking and later computation
            cross_velocity_squared_norm = (
                     (cross_velocities[0])*(cross_velocities[0])
                   + (cross_velocities[1])*(cross_velocities[1]))

            # If norm squared is not very small, we calculate time to 
            # CPA. Else we set it at previous value
            if (np.abs(cross_velocity_squared_norm) > TOL):

                # Calculating time to CPA per given formula
                time_to_CPA = -(
                    (distance_at_zero[0]*cross_velocities[0]
                    + distance_at_zero[1]*cross_velocities[1])
                    / cross_velocity_squared_norm)

                # Converting time to CPA to seconds
                time_to_CPA_sec = (time_to_CPA 
                                    *METERS_DEGREES_CONVERSION 
                                    / 0.514444)

                # Creating distance vector and computing the absolute value
                distvec = (distance_at_zero 
                        + time_to_CPA*cross_velocities)

                CPA_distance = np.sqrt(distvec[0]*distvec[0]
                                        + distvec[1]*distvec[1])

            elif ((np.abs(cross_velocity_squared_norm) < TOL)
                    &(index_ships != situation_start)):

                # Propagating past values
                time_to_CPA_sec = time_to_CPA_collector[index_ships-1]
                CPA_distance = CPA_dist_collector[index_ships-1]

            else:
                # Creating arbitrary values 
                time_to_CPA_sec = 9999999
                CPA_distance = 9999999

            # Storing the obtained values
            time_to_CPA_collector[index_ships] = time_to_CPA_sec
            CPA_dist_collector[index_ships] = CPA_distance

            # Calculating and storing distance to CPA
            d_to_CPA_ship1 = np.abs(speed1)*time_to_CPA_sec*0.514444
            d_to_CPA_ship2 = np.abs(speed2)*time_to_CPA_sec*0.514444

            dist_to_CPA1_collector[index_ships] = d_to_CPA_ship1
            dist_to_CPA2_collector[index_ships] = d_to_CPA_ship2

    # Correcting for the fact the CPA distance is calculated in terms of
    # degrees not meters
    CPA_dist_collector = CPA_dist_collector*METERS_DEGREES_CONVERSION

    return (time_to_CPA_collector, 
            CPA_dist_collector, 
            dist_to_CPA1_collector, 
            dist_to_CPA2_collector
            )


def statistics_aggregator(whole_table, 
                    timemergedcolname = "Merge_times", 
                    sitcolname = "Situations_1", 
                    idcolname1 = "ID_1",
                    idcolname2 = "ID_2",
                    SOGcolname1 = "SOG_1",
                    SOGcolname2 = "SOG_2",
                    COGcolname1 = "COG_1",
                    COGcolname2 = "COG_2",
                    LONcolname1 = "LON_1",
                    LONcolname2 = "LON_2",
                    LATcolname1 = "LAT_1",
                    LATcolname2 = "LAT_2",
                    dSOGcolname1 = "dSOG_1",
                    dSOGcolname2 = "dSOG_2",
                    dCOGcolname1 = "dCOG_1",
                    dCOGcolname2 = "dCOG_2",
                    tCOAcolname = "Time_to_CPA",
                    CPAdistcolname = "CPA_distance",
                    disttoCPAname1 = "Distance_to_CPA_1",
                    disttoCPAname2 = "Distance_to_CPA_2",
                    catcolname1 = "Category_1",
                    catcolname2 = "Category_2"
                    ):
    """
    Function to produce a table of summary statistics for each situation 
    collected from a table of merged, synchronized situation records.

    Inout: Synchronized table of situation records
    Output: Summary statistics table
    """

    # Sorting values in order for further computational choices to work.
    # Resetting index in order to make sure indexing works as it should
    whole_table = whole_table.sort_values([sitcolname,timemergedcolname])
    whole_table = whole_table.reset_index(drop=True)

    # Calling workhorse function
    (sitnum, recnum, t1_inds ,tF_inds ,t1s, tFs, tMs,
            delTimes, meanSOG1, meanSOG2, meanCOG1, meanCOG2,deltaCOG,
            COLREGs, dist_yields1, dist_yields2, pass_dist, app_speed, 
            tot_SOG1_change, tot_SOG2_change, tot_COG1_change, tot_COG2_change,
            max_dSOG1, max_dSOG2, max_dCOG1, max_dCOG2, 
            COG_maneuver1, COG_maneuver2, SOG_maneuver1, SOG_maneuver2, 
            ID1s, ID2s, Cat1, Cat2) = (
            _aggregator_workhorse(
        np.asarray(whole_table[idcolname1]),
        np.asarray(whole_table[idcolname2]),
        np.asarray(whole_table[LATcolname1]),
        np.asarray(whole_table[LATcolname2]),
        np.asarray(whole_table[LONcolname1]),
        np.asarray(whole_table[LONcolname2]),
        np.asarray(whole_table[timemergedcolname]),
        np.asarray(whole_table[SOGcolname1]),
        np.asarray(whole_table[SOGcolname2]),
        np.asarray(whole_table[COGcolname1]),
        np.asarray(whole_table[COGcolname2]),
        np.asarray(whole_table[sitcolname]),
        np.asarray(whole_table[disttoCPAname1]),
        np.asarray(whole_table[disttoCPAname2]),
        np.asarray(whole_table[dSOGcolname1]),
        np.asarray(whole_table[dSOGcolname2]),
        np.asarray(whole_table[dCOGcolname1]),
        np.asarray(whole_table[dCOGcolname2]),
        np.asarray(whole_table[tCOAcolname]),
        np.asarray(whole_table[CPAdistcolname]),
        np.asarray(whole_table[catcolname1]),
        np.asarray(whole_table[catcolname2])
        ))

    # Creating new output table (optional columns included in comments)
    outtable = pd.DataFrame({
            'Situation':sitnum, 'Records':recnum, 
            'Yield_start_record': t1_inds ,'Yield_finish_record': tF_inds,
            'ID_1': ID1s, 'ID_2': ID2s,"Category_1":Cat1, "Category_2":Cat2,
            'Start':t1s,'Finish' :tFs,'Midpoint':tMs,
            'Time_spent': delTimes, 'Seconds_spent':delTimes,
            #'Mean_SOG_1': meanSOG1,'Mean_SOG_2': meanSOG2,
            #'Mean_COG_1': meanCOG1,'Mean_COG_2': meanCOG2,
            'COLREG': COLREGs,'Crossing_angle': deltaCOG,
            'Distance_at_yield_1': dist_yields1,
            'Distance_at_yield_2': dist_yields2,
            'Passing_distance': pass_dist,
            'Approach_speed': app_speed, 
            'Max_speed_change_1': tot_SOG1_change,
            'Max_speed_change_2': tot_SOG2_change,
            'Max_course_change_1': tot_COG1_change,
            'Max_course_change_2': tot_COG2_change,
            #'Max_dSOG_1': max_dSOG1,'Max_dSOG_2': max_dSOG2,
            #'Max_dCOG_1': max_dCOG1,'Max_dCOG_2': max_dCOG2, 
            'Course_maneuver_1': COG_maneuver1,
            'Course_maneuver_2':COG_maneuver2,
            'Speed_maneuver_1': SOG_maneuver1,
            'Speed_maneuver_2': SOG_maneuver2},
            columns=['Situation','Records',
                    'Yield_start_record', 'Yield_finish_record',
                    'Start','Midpoint','Finish' ,
                    'Time_spent','Seconds_spent',
                    'ID_1','ID_2','Category_1','Category_2',
                    #'Mean_SOG_1', 'Mean_SOG_2',
                    #'Mean_COG_1', 'Mean_COG_2', 
                    'Crossing_angle', 'COLREG', 
                    'Distance_at_yield_1', 
                    'Distance_at_yield_2', 
                    'Passing_distance', 
                    'Approach_speed',  
                    'Max_speed_change_1', 
                    'Max_speed_change_2', 
                    'Max_course_change_1', 
                    'Max_course_change_2', 
                    #'Max_dSOG_1', 'Max_dSOG_2', 
                    #'Max_dCOG_1', 'Max_dCOG_2',  
                    'Course_maneuver_1', 'Course_maneuver_2',
                    'Speed_maneuver_1','Speed_maneuver_2'] 
            )

    # Selecting non-excluded records (choice of ID and Start as selectors is 
    # somewhat arbitrary)
    outtable = outtable[outtable.Start != 0]
    outtable = outtable[outtable.ID_1 != 0]
    
    # Converting time point values in integers (seconds) to datetime objects
    outtable.Start = pd.to_datetime(
                        outtable.Start*10**9).dt.floor('s')
    
    outtable.Midpoint = pd.to_datetime(
                        outtable.Midpoint*10**9).dt.floor('s')
    
    outtable.Finish = pd.to_datetime(
                        outtable.Finish*10**9).dt.floor('s')
    
    # Making a time delta value (second value is kept, see above)
    outtable.Time_spent = pd.to_timedelta(
                        outtable.Time_spent, unit="sec").dt.floor('s')

    # Resetting index 
    outtable = outtable.reset_index(drop=True)

    return outtable

def _aggregator_workhorse(ID1s, ID2s, LAT1s, LAT2s, LON1s, LON2s, Mergetimes, 
        SOG1s, SOG2s, COG1s, COG2s, Situations, dist2CPA1s, dist2CPA2s,
        dSOG1s, dSOG2s, dCOG1s, dCOG2s, tCPAs, CPAdists, Category1, Category2):
    """
    Workhorse function for calculating summary statistics by iterating through 
    each situation.
    """


    # Creating tolerance for small zero comparisons and creating
    # Conversion constant between degrees and meters
    TOL = 1E-7
    METERS_DEGREES_CONVERSION = 1.1132e+5

    # Asserting iteration values for situation loop
    situation_nums = np.unique(Situations)
    first_situation = np.uint32(situation_nums[0])
    number_of_sits = np.uint32(np.max(Situations) + 1)

    # Creating collection arrays for all relevant statistics (NB: Not all
    # arrays are output, but are kept in order for this to be a
    # future option)
    situation_numbers = np.array(range(0,number_of_sits),dtype=np.uint32)
    records_num_collector = np.zeros(number_of_sits,dtype = np.uint32)

    t1_ind_collector = np.zeros(number_of_sits,dtype = np.uint32)
    tF_ind_collector = np.zeros(number_of_sits,dtype = np.uint32)

    t1_collector = np.zeros(number_of_sits,dtype = np.float64)
    tF_collector = np.zeros(number_of_sits,dtype = np.float64)
    tM_collector = np.zeros(number_of_sits,dtype = np.float64)
    deltaTime_collector = np.zeros(number_of_sits,dtype = np.float64)

    mean_SOG1_collector = np.zeros(number_of_sits,dtype = np.float64)
    mean_SOG2_collector = np.zeros(number_of_sits,dtype = np.float64)
    mean_COG1_collector = np.zeros(number_of_sits,dtype = np.float64)
    mean_COG2_collector = np.zeros(number_of_sits,dtype = np.float64)

    delta_COG_collector = np.zeros(number_of_sits,dtype = np.float64)
    COLREG_collector = np.zeros(number_of_sits,dtype = np.uint32)

    distance_yield1_collector = np.zeros(number_of_sits,dtype = np.float64)
    distance_yield2_collector = np.zeros(number_of_sits,dtype = np.float64)
    passing_dist_collector = np.zeros(number_of_sits,dtype = np.float64)

    approach_speed_collector = np.zeros(number_of_sits,dtype = np.float64)

    tot_SOG1_change_collector = np.zeros(number_of_sits,dtype = np.float64)
    tot_SOG2_change_collector = np.zeros(number_of_sits,dtype = np.float64)
    tot_COG1_change_collector = np.zeros(number_of_sits,dtype = np.float64)
    tot_COG2_change_collector = np.zeros(number_of_sits,dtype = np.float64)

    max_dSOG1_collector = np.zeros(number_of_sits,dtype = np.float64)
    max_dSOG2_collector = np.zeros(number_of_sits,dtype = np.float64)
    max_dCOG1_collector = np.zeros(number_of_sits,dtype = np.float64)
    max_dCOG2_collector = np.zeros(number_of_sits,dtype = np.float64)

    COG_maneuver1_collector = np.zeros(number_of_sits,dtype = np.bool_)
    COG_maneuver2_collector = np.zeros(number_of_sits,dtype = np.bool_)
    SOG_maneuver1_collector = np.zeros(number_of_sits,dtype = np.bool_)
    SOG_maneuver2_collector = np.zeros(number_of_sits,dtype = np.bool_)

    ID_collector1 = np.zeros(number_of_sits,dtype = np.uint32)
    ID_collector2 = np.zeros(number_of_sits,dtype = np.uint32)

    Cat1_collector1 = np.zeros(number_of_sits,dtype = np.uint32)
    Cat1_collector2 = np.zeros(number_of_sits,dtype = np.uint32)

    # Looping through all situations
    for situation in range(first_situation,first_situation + number_of_sits):

        # Asserting iteration values for in-situation loop
        situation_filter = Situations == situation
        records_num = np.sum(situation_filter)

        # Storing number of records
        records_num_collector[situation] = records_num

        # Exiting "empty" situations
        if records_num == 0:
            pass
        else:
            
            # Retrieving CPA_distances and Times to CPA for all
            # records
            CPAdist_sit = CPAdists[situation_filter]
            tCPA_sit = tCPAs[situation_filter]

            # Quitting uninteresting situations with too high CPA_dist and
            # where Time to CPA is negative throughout the situation interval
            if (tCPA_sit<=0).all():
                pass

            elif (CPAdist_sit > 50).all():
                pass

            else:
                ### Finding time intervals ###

                # Setting off a separate merged time interval to work on
                Mergetime_sit = Mergetimes[situation_filter]

                # Option one: CPA_dist never goes below 10. Then we select the
                # instance where tCPA is positive and CPA_distance is smallest. 
                # If this does not exist, we exit and leave the situation (not
                # implemented)
                if (CPAdist_sit >= 10).all():
                    if ((CPAdist_sit >= 10)&(tCPA_sit>0)).any():
                        t1_ind = np.argmin(CPAdist_sit)

                    else:
                        msg = "You have not taken account of this situation"
                        raise NotImplementedError(msg)

                else:
                    # If there are instances with CPA distance below 10, we 
                    # proceed with finding the first instance where CPA distance 
                    # goes below 10 meters, and Time to CPA is positive.
                    if ((CPAdist_sit < 10)&(tCPA_sit>0)).any():

                        t1_ind_filter = np.logical_and(
                                            CPAdist_sit < 10,
                                            tCPA_sit>0)

                        t1_ind = np.where(t1_ind_filter)
                        t1_ind = t1_ind[0][0]

                    # If CPA distance goes below 10, but does not for a record 
                    # where Time to CPA is positive, we find the smallest
                    # CPA_distance value where Time to CPA is positive, as
                    # above. If this does not exist, we exit and leave the 
                    # situation (not implemented)
                    elif ((CPAdist_sit >= 10)&(tCPA_sit>0)).any():
                        t1_ind = np.argmin(CPAdist_sit)

                    else:
                        msg = "You have not taken account of this situation"
                        raise NotImplementedError(msg)

                t1 = Mergetime_sit[t1_ind]
                
                # We then proceed to find the time point where the vessel
                # passes out of the near-collision situation. We find this as
                # the first record after t1 where Time to CPA turns negative.
                # If this does not exist, we choose the last record in the set
                # as our exit point (implemented in reverse order)
                if not np.logical_and(tCPA_sit < 0, Mergetime_sit > t1).any():
                    tF_ind = records_num - 1

                else:                
                    tF_ind_filter = np.logical_and(
                                            tCPA_sit < 0, 
                                            Mergetime_sit > t1)

                    tF_ind = np.where(tF_ind_filter)
                    tF_ind = tF_ind[0][0]

                # Asserting finish time using index. Calculating midpoint time
                # as simple aritmetric mean
                tF = Mergetime_sit[tF_ind]
                tM = (t1 + tF)/2

                # Checking if there is more than 2 records between start point
                # and midpoint. Exiting if not.
                if (np.sum(np.logical_and(Mergetime_sit >= t1 , 
                                        Mergetime_sit <= tM)) <= 2):
                    pass

                else:

                    # Collecting IDs and categories. The main function uses this
                    # assertion to filter for exited situations later
                    ID_collector1[situation] = np.unique(
                                    ID1s[situation_filter])[0]
                    
                    ID_collector2[situation] = np.unique(
                                    ID2s[situation_filter])[0]
                    
                    Cat1_collector1[situation] = np.unique(
                                    Category1[situation_filter])[0]
                    
                    Cat1_collector2[situation] = np.unique(
                                    Category2[situation_filter])[0]

                    # Finding the time spent in the yielding maneuver
                    deltaTime = tF - t1

                    # Storing the computed time values
                    t1_collector[situation] = t1
                    tF_collector[situation] = tF
                    tM_collector[situation] = tM 
                    deltaTime_collector[situation] = deltaTime

                    # Storing indices for robustness analysis of final data
                    t1_ind_collector[situation] = t1_ind 
                    tF_ind_collector[situation] = tF_ind

                    ### Finding distances ###

                    # Collecting Distances to CPA at point of yield for both
                    # ships (calculated earlier), and passing distance as
                    # CPA_dist at finishing time of yield maneuver
                    distance_yield1_collector[situation] = (
                                        dist2CPA1s[situation_filter][t1_ind])

                    distance_yield2_collector[situation] = (
                                        dist2CPA2s[situation_filter][t1_ind])

                    passing_dist_collector[situation] = CPAdist_sit[tF_ind]

                    ### Categorizing situations in COLREGs ###

                    # Collecting COGs and SOGs for situation as separate array
                    COG1s_sit = COG1s[situation_filter]
                    COG2s_sit = COG2s[situation_filter]
                    SOG1s_sit = SOG1s[situation_filter]
                    SOG2s_sit = SOG2s[situation_filter]

                    # Computing mean SOGs and COGs as mean of COGs and SOGs 
                    # between start of yield and midpoint of yield maneuver
                    mean_COG1_firsthalf = np.mean(COG1s_sit[
                        np.logical_and(Mergetime_sit >= t1 , 
                                        Mergetime_sit <= tM)])
                    
                    mean_COG2_firsthalf = np.mean(COG2s_sit[
                        np.logical_and(Mergetime_sit >= t1 , 
                                        Mergetime_sit <= tM)])
                    
                    mean_SOG1_firsthalf = np.mean(SOG1s_sit[
                        np.logical_and(Mergetime_sit >= t1 , 
                                        Mergetime_sit <= tM)])
                    
                    mean_SOG2_firsthalf = np.mean(SOG2s_sit[
                        np.logical_and(Mergetime_sit >= t1 , 
                                        Mergetime_sit <= tM)])

                    # Saving mean_COGs and mean_SOGs for output if needed
                    # (not currently output)
                    mean_SOG1_collector[situation] = mean_SOG1_firsthalf
                    mean_SOG2_collector[situation] = mean_SOG2_firsthalf
                    mean_COG1_collector[situation] = mean_COG1_firsthalf
                    mean_COG2_collector[situation] = mean_COG2_firsthalf


                    # Finding the mean COG from before entering into the
                    # yielding maneuver for both vessels
                    mean_COG1_pre = np.mean(COG1s_sit[(Mergetime_sit <= t1)
                                                &(Mergetime_sit > (t1 - 20))])

                    mean_COG2_pre = np.mean(COG2s_sit[(Mergetime_sit <= t1)
                                                &(Mergetime_sit > (t1 - 20))])

                    # Finding crossing angles by using the collected mean
                    # COG and SOGs and collecting this
                    delta_COG = (mean_COG1_pre - mean_COG2_pre) % 360
                    delta_COG_collector[situation] = delta_COG

                    # Categorizing into COLREG situation (1 is overtaking, 
                    # 2 is crossing and 3 is head-to-head) by using 
                    # the mean COG_pre found previously)
                    if np.cos(delta_COG*np.pi/180) > np.cos(45*np.pi/180):
                        COLREG_collector[situation] = 1 

                    elif np.cos(delta_COG*np.pi/180)< - np.cos(15*np.pi/180):
                        COLREG_collector[situation] = 3 

                    else:
                        COLREG_collector[situation] = 2 

                    ### Calculating approach speed ###

                    # Calculating approach velocities in LON and LAT direction
                    approach_speed_x = (
                            mean_SOG1_firsthalf*np.sin(
                                mean_COG1_firsthalf*np.pi/180) 
                            - mean_SOG2_firsthalf*np.sin(
                                mean_COG2_firsthalf*np.pi/180))

                    approach_speed_y = (
                            mean_SOG1_firsthalf*np.cos(
                                mean_COG1_firsthalf*np.pi/180) 
                            - mean_SOG2_firsthalf*np.cos(
                                mean_COG2_firsthalf*np.pi/180))

                    # Computing and collecting approach speed 
                    approach_speed_collector[situation] = np.sqrt(
                                approach_speed_x*approach_speed_x 
                                + approach_speed_y*approach_speed_y)

                    ### Calculating total change in COG and SOG during yield ###

                    # Retrieving relevant COGs and SOGs for the evasion maneuver
                    # interval
                    COG1_int = COG1s_sit[np.logical_and(
                                    Mergetime_sit >= t1 , 
                                    Mergetime_sit <= tF)]
                    
                    COG2_int = COG2s_sit[np.logical_and(
                                    Mergetime_sit >= t1 , 
                                    Mergetime_sit <= tF)]
                    
                    SOG1_int = SOG1s_sit[np.logical_and(
                                    Mergetime_sit >= t1 , 
                                    Mergetime_sit <= tF)]
                    
                    SOG2_int = SOG2s_sit[np.logical_and(
                                    Mergetime_sit >= t1 , 
                                    Mergetime_sit <= tF)]

                    # Computing maximum change in COG1 and COG2
                    tot_COG1_change = np.max(COG1_int - np.min(COG1_int))
                    tot_COG2_change = np.max(COG2_int - np.min(COG2_int))

                    # Correcting the COG change if a high number results from
                    # passing the 0/360 mark and giving a value above 180
                    # for both COGs
                    if tot_COG1_change > 180:
                        tot_COG1_change_collector[situation] = (
                                                    360 - tot_COG1_change)
                    else:
                        tot_COG1_change_collector[situation] = tot_COG1_change

                    if tot_COG2_change > 180:
                        tot_COG2_change_collector[situation] = (
                                                    360 - tot_COG2_change)
                    else:
                        tot_COG2_change_collector[situation] = tot_COG2_change

                        # Saving COG change values
                    tot_SOG1_change_collector[situation] = np.max(
                                                    SOG1_int - np.min(SOG1_int))
                    tot_SOG2_change_collector[situation] = np.max(
                                                    SOG2_int - np.min(SOG2_int))

                    ### Finding the type of manouver ###

                    # Retrieving dSOGs and dCOGs for the situation
                    dCOG1s_sit = dCOG1s[situation_filter]
                    dCOG2s_sit = dCOG2s[situation_filter]
                    dSOG1s_sit = dSOG1s[situation_filter]
                    dSOG2s_sit = dSOG2s[situation_filter]

                    # Finding largest differentiated COG and SOG before yield 
                    # maneuver
                    max_dCOG1_pre = np.max(dCOG1s_sit[Mergetime_sit <= t1])
                    max_dCOG2_pre = np.max(dCOG2s_sit[Mergetime_sit <= t1])
                    max_dSOG1_pre = np.max(dSOG1s_sit[Mergetime_sit <= t1])
                    max_dSOG2_pre = np.max(dSOG2s_sit[Mergetime_sit <= t1])

                    # Saving the largest dSOGs and dCOGs
                    max_dCOG1_collector[situation] = max_dCOG1_pre
                    max_dCOG2_collector[situation] = max_dCOG2_pre
                    max_dSOG1_collector[situation] = max_dSOG1_pre
                    max_dSOG2_collector[situation] = max_dSOG2_pre

                    # Retrieving dSOGs and dCOGs for the yielding interval
                    dCOG1s_int = dCOG1s_sit[np.logical_and(
                                            Mergetime_sit >= t1 , 
                                            Mergetime_sit <= tF)]
                    
                    dCOG2s_int = dCOG2s_sit[np.logical_and(
                                            Mergetime_sit >= t1 , 
                                            Mergetime_sit <= tF)]
                    
                    dSOG1s_int = dSOG1s_sit[np.logical_and(
                                            Mergetime_sit >= t1 , 
                                            Mergetime_sit <= tF)]
                    
                    dSOG2s_int = dSOG2s_sit[np.logical_and(
                                            Mergetime_sit >= t1 , 
                                            Mergetime_sit <= tF)]

                    # Stamping ship 1 and 2 with "True" if they had a maximum
                    # dSOG and dCOG exceeding the value had before the
                    # yielding maneuver started
                    if len(dCOG1s_int) > 0:
                        COG_maneuver1_collector[situation] = np.any(
                                            dCOG1s_int > max_dCOG1_pre)

                    if len(dCOG2s_int) > 0:
                        COG_maneuver2_collector[situation] = np.any(
                                            dCOG2s_int > max_dCOG2_pre)

                    if len(dSOG1s_int) > 0:
                        SOG_maneuver1_collector[situation] = np.any(
                                            np.abs(dSOG1s_int) > max_dSOG1_pre)

                    if len(dSOG2s_int) > 0:
                        SOG_maneuver2_collector[situation] = np.any(
                                            np.abs(dSOG2s_int) > max_dSOG2_pre)
            

    # Returning all records (not neccessarily storing all in final output of
    # main function)
    return (situation_numbers,
            records_num_collector,
            t1_ind_collector,
            tF_ind_collector,
            t1_collector,
            tF_collector,
            tM_collector,
            deltaTime_collector,
            mean_SOG1_collector,
            mean_SOG2_collector,
            mean_COG1_collector,
            mean_COG2_collector,
            delta_COG_collector,
            COLREG_collector,
            distance_yield1_collector,
            distance_yield2_collector,
            passing_dist_collector,
            approach_speed_collector,
            tot_SOG1_change_collector,
            tot_SOG2_change_collector,
            tot_COG1_change_collector,
            tot_COG2_change_collector,
            max_dSOG1_collector,
            max_dSOG2_collector,
            max_dCOG1_collector,
            max_dCOG2_collector,
            COG_maneuver1_collector,
            COG_maneuver2_collector,
            SOG_maneuver1_collector,
            SOG_maneuver2_collector,
            ID_collector1,
            ID_collector2,
            Cat1_collector1,
            Cat1_collector2
            )

if __name__ == "__main__":
    pass

