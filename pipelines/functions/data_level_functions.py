# datalevel_functions.py
"""
Description    :    Processing of data into various levels
Author         :    Karla Onate Melecio (KOM) <kgom.astro at gmail.com>
Date Created   :    2023/08/10
Date Modified  :    2024/05/09
Version        :    5.0
Python Version :    3.10.8   
"""

import sys
import struct
import pathlib
import os
import re
import logging
import warnings

import pandas as pd
import numpy as np
import uproot


########### Raw Data Level #######################################

# The raw data is first pre-processed with the convert_raw_data_format function which process the raw .ldat file from the Petsys system into a pandas DataFrame containing the unprocessed raw data. convert_raw_data_format utilizes a conversion function called convert_binary_to_dictionary. The convert_binary_to_dictionary function processes a binary file containing petsys "singles" data (all triggers) by reading it in chunks of 16 bytes, each representing a timestamp (ps time tag relative to system start time), energy value (digitized measure of charge), and 'channelID' (detector identification number). It then unpacks these values and stores them in lists. Once all data is read, it constructs a dictionary with keys 'time', 'channelID', and 'energy', where each key corresponds to the respective list of values. The convert_raw_data_format function, offers the flexibility to specify the output format, including Feather, compressed CSV, or compressed Pickle files. This step converts binary data into a structured format and saves it in various file formats specified by the user for further analysis or storage.

def convert_binary_to_dictionary(path_to_file):
    """
    Converts processed petsys-singles file to a dictionary of raw triggers.

    Args:
        path_to_file (str): Path to the binary file to be processed.

    Returns:
        dict: Dictionary containing time, channelID, and energy lists.
    """
    with open(path_to_file, "rb") as f:
        nbytes = 16 

        # Initializes lists for trigger info
        time_col = []
        chan_col = []
        energy_col = []

        # Unpacks the binary file and appends each value to the corresponding list
        while True:
            buffer = f.read(nbytes)
            if len(buffer) == 0:
                break
            elif len(buffer) == nbytes:
                time, energy, channelID = struct.unpack("<qfi", buffer)
                time_col.append(time)
                energy_col.append(energy)
                chan_col.append(channelID)
            else:
                print("Warning: Incomplete data detected")

    raw_triggers = {'time': time_col, 'channelID': chan_col, 'energy': energy_col}
    return raw_triggers


def convert_raw_data(run_id, directory='../data/raw/', output_format=None, input_format='.ldat'):
    """
    Converts .ldat or .root file to another specified format.

    Arguments:
        run_id (str): Filename of the .ldat file.
        directory (str): Path to the directory containing the .ldat file or .root.
        output_format (str): Desired output format ('.feather', 'csv.zip', '.pkl.gzip').

    Returns:
        pandas.DataFrame: DataFrame containing raw data.
    """
    if input_format == '.root':
        raw_df = UprootIt(os.path.join(directory, f'{run_id}.root'), outType="pd", cols=['time', 'channelID', 'energy'], cuts=None, alias=None, filters=[None, None, None])
    elif input_format == '.ldat':
        # Takes in raw .ldat file
        raw_df = pd.DataFrame(convert_binary_to_dictionary(os.path.join(directory, f'{run_id}.ldat')))

    # Produces one of three specified file formats containing raw data and adds to raw folder otherwise returns only pandas.DataFrame
    if output_format == '.feather':
        raw_df.to_feather(os.path.join(directory, f'{run_id}.feather'))
        print('Feather file created')
    elif output_format == 'csv.zip':
        raw_df.to_csv(os.path.join(directory, f'{run_id}.csv.zip'), index=False, compression='zip')
        print('CSV.zip file created')
    elif output_format == '.pkl.gzip':
        raw_df.to_pickle(os.path.join(directory, f'{run_id}.pkl.gzip'), compression='gzip')
        print('Pickle gzip file created')
    else:
        print('Supported file format not provided or invalid format specified. Returning pandas.DataFrame; no data file created')

    return raw_df

# The process_raw_data_directory function iterates through the .ldat files in a specified directory. For each file, it extracts the run ID, invokes the convert_raw_data function to process the binary file into a structured DataFrame, and saves it in the desired output format.

def process_raw_data_directory(directory='../data/raw/', output_format='.feather', input_format='.ldat'):
    """
        Processes raw data files from the specified directory. For each file, it extracts the run ID, invokes the convert_raw_data function to process the binary file into a structured DataFrame, and saves it in the desired output format.

        Parameters:
        directory (str): Path to the directory containing the .ldat file.
        output_format (str): Desired output format ('.feather', 'csv.zip', '.pkl.gzip').

        Returns:
        None
    """
    # gets new datafiles with picosecond resolution
    data = pathlib.Path(directory)

    # creates a list of the paths to the binary files in the data folder
    files = list(data.glob("*.ldat"))

    for i in files:
        s = str.split(str(i), sep='raw\\')[1]
        run_id = re.search('(.*).ldat', s).group(1)

        print(run_id)

        convert_raw_data(run_id, directory=directory, output_format=output_format, input_format=input_format)

    return



############################################################################################################
# UprootIt: warnings, uproot, pandas
############################################################################################################
"""
UprootIt takes in a root file and converts it to the desired output format (np, pd, ak)
UprootIt(fInfo, outType, cols, cuts, alias, filters)
    fInfo : Expects file name and path as a string
    outType:filters are select arguments for uproot TTree.arrays
    outType: ("pd", "np", "ak") this is the library argument for the uproot TTree.arrays function. Default value is "pd"
"""
def UprootIt(fInfo, outType="pd", cols=None, cuts=None, alias=None, filters=[None, None, None]):

    tree = uproot.open(fInfo + ':data')

    # suppresses copy induced performance warnings temporarily
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning) 
    
    nam = filters[0]
    typ = filters[1]
    bch = filters[2]

    df = tree.arrays(expressions=cols, cut=cuts, filter_name=nam, filter_typename=typ, filter_branch=bch, aliases=alias, entry_start=None, entry_stop=None, library=outType)
    
    warnings.simplefilter(action='default', category=pd.errors.PerformanceWarning)

    return df

##################################################################


def convert_raw_adc_to_adc(raw_adc_df):
    '''
    Converts raw adc data to ASIC corrected (ADC) values.
    
    Arguments:
        raw_adc_df (pandas.DataFrame): DataFrame containing raw adc data.
        
    Returns:
        pandas.DataFrame: DataFrame containing ADC values.
    '''
    # Coefficients for ADC conversion
    a = 8.00000
    b = 1.04676
    c = 1.02734
    d = 0.31909

    # Calculate ADC values based on the formula for each energy value
    Q = raw_adc_df.energy.copy()
    raw_adc_df['adc'] =  ((a)*((b)**(Q**(c)))) + (((d)*Q)-(a))

    # Display the number of triggers before filtering
    logging.info(f'Unfiltered input contains {len(Q)} triggers')

    # Remove the 'energy' column and return the DataFrame with ADC values
    adc_df = raw_adc_df.drop(columns=['energy']).copy()
    return adc_df


def remove_adc_glitches(adc_df, new_glitch_min=None, new_glitch_max=None):
    '''
    Removes channel glitches from ADC data.
    
    Arguments:
        adc_df (pandas.DataFrame): DataFrame containing ADC values.
        new_glitch_min (float or int): Value corresponding to a new minimum of channel interval to be removed
        new_glitch_max (float or int): Value corresponding to a new maximum of channel interval to be removed
        
    Returns:
        pandas.DataFrame: DataFrame with known ADC channel glitches removed.
    '''
    # Filtering ADC values based on predefined intervals
    adc_df = adc_df[(adc_df.adc <= 182.92) | (adc_df.adc >= 182.93)]
    adc_df = adc_df[(adc_df.adc <= 510.63) | (adc_df.adc >= 510.64)]
    adc_df = adc_df[adc_df.adc < 1000]

    # Remove glitches within specified intervals
    if (new_glitch_min is not None) and (new_glitch_max is not None):
        adc_df = adc_df[(adc_df.adc <= new_glitch_min) | (adc_df.adc >= new_glitch_max)]

    # Display the number of triggers after filtering
    logging.info(f'Filtered output contains {len(adc_df)} triggers')
    return adc_df    


def L0_data_preprocess(input_data, input_data_format, output_format, run_id=None):
    '''
    Pre-processes raw data to L0 data level.

    Arguments:
        input_data (str or pandas.DataFrame): Input data, either file path, DataFrame, or .ldat.
        input_data_format (str): Type of input data ('.pkl.gzip','.csv.zip', '.feather', 'dataframe', '.root').
        output_format (str): Desired output format ('.feather', '.csv.zip', '.pkl.gzip').
        run_id (str, optional): Identifier for the run. Required if output_format is not empty.

    Returns:
        pandas.DataFrame: Cleaned DataFrame.
    '''
    # Load raw data based on the input data format
    if input_data_format == '.pkl.gzip':
        raw_adc_df = pd.read_pickle(input_data, compression='gzip')
    elif input_data_format == '.csv.zip':
        raw_adc_df = pd.read_csv(input_data, compression='zip')
    elif input_data_format == '.feather':
        raw_adc_df = pd.read_feather(input_data)
    elif input_data_format == 'dataframe':
        raw_adc_df = input_data
    elif input_data_format == '.root':
        raw_adc_df = UprootIt(os.path.join(directory, f'{run_id}.root'), outType="pd", cols=['time', 'channelID', 'energy'], cuts=None, alias=None, filters=[None, None, None])
    else:
        raise ValueError("Invalid input_data_format. Specify: ('.pkl.gzip','.csv.zip', '.feather', 'dataframe', '.root')")

    logging.info('raw_adc_df created')

    # Process data: convert ADC channels to ADC values
    if 'energy' in raw_adc_df.columns:
        adc_df = convert_raw_adc_to_adc(raw_adc_df)
        logging.info('raw ADC channels converted to ADC values')
    else:
        adc_df = raw_adc_df

    # Filter out known ADC channel glitches
    if 'adc' in adc_df.columns:
        clean_df = remove_adc_glitches(adc_df)
        clean_df = clean_df.reset_index(drop=True)
        logging.info('L0 data pre-processed')

    # Save processed data if run_id is provided
    if run_id is not None:
        output_path = f'../data/L0/L0_{run_id}'
        if output_format is not None:
            if os.path.exists(output_path + output_format):
                logging.warning('Output file already exists. Overwriting...')
        
        if output_format == '.feather':
            clean_df.to_feather(output_path + '.feather')
            logging.info('L0 .feather file created')
        elif output_format == '.csv.zip':
            clean_df.to_csv(output_path + '.csv.zip', index=False, compression='zip')
            logging.info('L0 ".csv.zip" file created')
        elif output_format == '.pkl.gzip':
            clean_df.to_pickle(output_path + '.pkl.gzip', compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})
            logging.info('L0 .pkl.gzip file created')
        elif output_format is not None:
            raise ValueError('Unsupported output file format')
    elif output_format is not None:
        raise ValueError('run_id is required for saving the output')

    return clean_df


def process_L0_data_from_raw(raw_directory='../data/raw/', L0_directory='../data/L0', input_format='.feather', output_format='.feather'):
    """
        Processes L0 data files from the raw directory, adds ancillary information (columnID, scintillatorType, PosZ, PosX, PosY),
        and saves the processed data into the L0 directory in specified format.

        Parameters:
        raw_directory (str): Path to the directory containing raw data files.
        L0_directory (str): Path to the directory where L0 data files will be saved.
        input_format (str): File extension of raw files.
        output_format (str): File extension of processed files.

        Returns:
        None
    """
    # Get the paths to all raw data files
    data = pathlib.Path(raw_directory)
    files = list(data.glob(f"*{input_format}"))

    # Process each raw data file
    for i in files:
        # Extract run ID from the file name
        rid = os.path.splitext(os.path.basename(i))[0]

        # Preprocess raw data to L0 level and save
        L0_data_preprocess(input_data=i, input_data_format=input_format, 
                           output_format=output_format, run_id=rid)

    return



######################################################################
# The process_L0_file function extracts L0 data file and returns the processed dataframes along with the run ID. The identify_events function processes events from the time-sorted data based on a specified time window, identifies singles events, and generates an event lookup DataFrame. The preprocess_L2_data function preprocesses L0 data files to extract singles events, saves the processed singles data, and event lookup information. It iterates over all L0 data files in a directory, processes each file, and saves the results to specified output directories.

def process_L0_file(file):
    """
    Process an L0 file to extract singles events before adding additional info.

    Parameters:
    file (pathlib.Path): Path to the L0 data file.

    Returns:
    Tuple[DataFrame, str]: Processed L0 data with singles events extracted and run ID.
    """
    # Extract run ID from the file name
    rid = re.search(r'L0_(.*).feather', file.name).group(1)
    logging.info(f"Processing run ID: {rid}")
    
    # Read the Feather file into a DataFrame
    data_0 = pd.read_feather(file)

    # Sort the DataFrame by time
    data_0.sort_values(by='time', inplace=True)
    
    # Reset the index
    data_0.reset_index(drop=True, inplace=True)
    
    return data_0, rid

def identify_events(data, window=4000):
    """
    Process events from time-sorted data based on a specified time window.

    Parameters:
    data (DataFrame): The input data containing time values.
    window (int): The time window in picoseconds for singles multiplicity (default is 10000).

    Returns:
    Tuple[DataFrame, DataFrame]: DataFrame with singles events and event lookup DataFrame.
    """
    # Initialize lists to store events and their characteristics
    event_entries = []
    first_hit = []
    multiplicity = []

    # Initialize the first event
    event_start = data['time'].iloc[0]
    current_event = []
    current_event_start_index = 0

    logging.info(f"Identifying events")
    # Loop through the time values in the data with a progress bar
    for i, hit_time in enumerate(data['time']):
        # Check if the time difference between the current hit and the start of the event exceeds the window
        if hit_time - event_start >= window:
            # End the current event and start a new one
            event_entries.append(current_event)
            multiplicity.append(len(current_event))
            first_hit.append(current_event_start_index)

            # Start a new event
            current_event = []
            event_start = hit_time
            current_event_start_index = i

        # Add the current hit to the current event
        current_event.append(i)

    # Append the last event after the loop
    event_entries.append(current_event)
    multiplicity.append(len(current_event))
    first_hit.append(current_event_start_index)

    # Create a new column 'event_id' in data for efficient filtering      
    data["event_id"] = data.index

    # Create event lookup DataFrame
    event_lookup_df = pd.DataFrame({
        "event_id": first_hit,
        "multiplicity": multiplicity,
        "hits": event_entries
    })

    logging.info(f"Preparing singles events")
    # Filter data based on multiplicity using the event_id
    # singles_data = data[data.event_id.isin(event_lookup_df[event_lookup_df.multiplicity == 1].event_id)]

    logging.info(f"Preparing event lookup table")
    # Expand the 'event_id' to cover all hits in each event
    event_lookup_df = event_lookup_df.explode('hits').reset_index(drop=True)

    return event_lookup_df

def preprocess_L2_data(directory, output_directory, lookup_directory, window=4000):
    """
    Preprocess L0 data to extract singles events and save the processed data and event lookup.

    Parameters:
    directory (str): Path to the directory containing L0 data files.
    output_directory (str): Path to the directory to save processed singles data.
    lookup_directory (str): Path to the directory to save event lookup data.

    Returns:
    None
    """
    files = list(pathlib.Path(directory).glob("*.feather"))

    if not files:
        raise FileNotFoundError(f"No feather files found in directory {directory}")

    for file in files:
        data_0, run_id = process_L0_file(file)
        event_lookup_df = identify_events(data_0, window=window)

        # Reset indexes
        # singles.reset_index(drop=True, inplace=True)
        event_lookup_df.reset_index(drop=True, inplace=True)

        logging.info(f"Saving processed data for run ID: {run_id}")
        # # Save processed data
        # singles.to_feather(f"{output_directory}/L2A1_{run_id}.feather")
        event_lookup_df.to_feather(f"{lookup_directory}/events_{window}ps_window_{run_id}.feather")








########################### General ##############################

def measure_execution_time(func, *args, **kwargs):
    '''
    Measures the execution time of a function.

    Arguments:
        func (function): The function to be executed.
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Returns:
        tuple: Result of the function and execution time.
    '''
    start_time = time.time()
    result = func(*args, **kwargs)
    end_time = time.time()
    execution_time = end_time - start_time
    return result, execution_time

##################################################################


# def add_ancillaries_to_data(raw_triggers, identifying_data='ancillary_data/identifying_data.csv', drop_cols=['biasID', 'pinID']):
#     """
#         add_ancillaries_to_data _summary_

#     Parameters
#     ----------
#         raw_triggers : dictionary or dataframe from petsys containing the time, energy, and detector information
#         identifying_data : str, optional
#             ancillary lookup table with 'channelID' key, by default 'ancillary_data/identifying_data.csv'
#         drop_cols : list, optional
#             list of columns to drop if any, by default ['biasID', 'pinID']

#     Returns
#     -------
#         pandas dataframe 
#             pandas dataframe containing input data and ancillary data
#     """

#     trigger_df = pd.DataFrame(raw_triggers.copy())
#     identifying_df = pd.read_csv(identifying_data)

#     coded_data = pd.merge(trigger_df, identifying_df, how='left', on='channelID')
#     coded_data = coded_data.drop(columns=drop_cols)

#     return coded_data


# def convert_raw_data(run_id, output_format=None):
#     '''
#     Converts .ldat (.ldat) file to another specified format.
    
#     Arguments:
#         run_id (str): Filename of the .ldat file.
#         output_format (str): Desired output format ('.feather', 'csv.zip', '.pkl.gzip').
        
#     Returns:
#         pandas.DataFrame: DataFrame containing raw data.
#     '''
#     # Takes in raw .ldat file
#     raw_df = pd.DataFrame(convert_binary_to_dictionary('../data/raw/'+ run_id + '.ldat'))

#     # Produces one of three specified file formats containing raw data and adds to raw folder otherwise returns only pandas.DataFrame
#     if output_format == '.feather':
#         raw_df.to_feather('../data/raw/'+run_id+'.feather')
#         print('feather file created')
#     elif output_format == 'csv.zip':
#         raw_df.to_csv('../data/raw/'+run_id+'.csv.zip', index=False, compression='zip')
#         print('.csv.zip file created')
#     elif output_format == '.pkl.gzip':
#         raw_df.to_pickle('../data/raw/'+ run_id +'.pkl.gzip', compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})
#         print('.pkl.gzip file created')
#     elif output_format == None:
#         print('Supported file format not provided. Returning pandas.DataFrame; no data file created')

#     return raw_df

# #####

# def convert_raw_adc_to_adc(raw_adc_df):
#     '''
#     Converts raw adc data to ASIC corrected (ADC) values.
    
#     Arguments:
#         raw_adc_df (pandas.DataFrame): DataFrame containing raw adc data.
        
#     Returns:
#         pandas.DataFrame: DataFrame containing ADC values.
#     '''
#     a = 8.00000
#     b = 1.04676
#     c = 1.02734
#     d = 0.31909

#     Q = raw_adc_df.energy.copy()
#     raw_adc_df['adc'] =  ((a)*((b)**(Q**(c)))) + (((d)*Q)-(a))

#     print(f'Unfiltered input contains {len(Q)} triggers')
#     adc_df = raw_adc_df.drop(columns=['energy']).copy()

#     return adc_df


# def remove_adc_glitches(adc_df, new_glitch_min=None, new_glitch_max=None):
#     '''
#     Removes channel glitches from ADC data.
    
#     Arguments:
#         adc_df (pandas.DataFrame): DataFrame containing ADC values.
#         new_glitch_min (float or int): Value corresponding to a new minimum of channel interval to be removed
#         new_glitch_max (float or int): Value corresponding to a new maximum of channel interval to be removed
        
#     Returns:
#         pandas.DataFrame: DataFrame with known ADC channel glitches removed.
#     '''
#     # Filtering ADC values
#     adc_df = adc_df[(adc_df.adc <= 182.92) | (adc_df.adc >= 182.93)]
#     adc_df = adc_df[(adc_df.adc <= 510.63) | (adc_df.adc >= 510.64)]
#     adc_df = adc_df[adc_df.adc < 1000]

#     if (new_glitch_min is not None) and (new_glitch_max is not None):
#         adc_df = adc_df[(adc_df.adc <= new_glitch_min) | (adc_df.adc >= new_glitch_max)]

    
#     print(f'Filtered output contains {len(adc_df)} triggers')
#     return adc_df


# def L0_data_preprocess(input_data_format, input_data, output_format, run_id=None):
#     '''
#     Preprocesses raw data to L0 data level.

#     Arguments:
#         input_data (str or pandas.DataFrame): Input data, either file path, DataFrame, or .ldat.
#         input_data_format (str): Type of input data ('.pkl.gzip','.csv.zip', '.feather', 'dataframe').
#         output_format (str): Desired output format ('.feather', '.csv.zip', '.pkl.gzip').
#         run_id (str, optional): Identifier for the run. Required if output_format is not empty.

#     Returns:
#         pandas.DataFrame: Cleaned DataFrame.
#     '''
#     # Load raw data
#     if input_data_format == '.pkl.gzip':
#         raw_adc_df = pd.read_pickle(input_data, compression='gzip')
#     elif input_data_format == '.csv.zip':
#         raw_adc_df = pd.read_csv(input_data, compression='zip')
#     elif input_data_format == '.feather':
#         raw_adc_df = pd.read_feather(input_data)
#     elif input_data_format == 'dataframe':
#         raw_adc_df = input_data

#     else:
#         raise ValueError("Invalid input_data_format. Specify: ('.pkl.gzip','.csv.zip', '.feather', 'dataframe')")

#     logging.info('raw_adc_df created')

#     # Process data 
#     if 'energy' in raw_adc_df.columns:
#         adc_df = convert_raw_adc_to_adc(raw_adc_df)
#         logging.info('raw ADC channels converted to adc values')
#     else:
#         adc_df = raw_adc_df

#     if 'adc' in adc_df.columns:
#         clean_df = remove_adc_glitches(adc_df)
#         if 'idx' in clean_df.columns:
#             clean_df = clean_df.drop(columns=['idx'])
#         elif 'index' in clean_df.columns:
#             clean_df = clean_df.drop(columns=['index'])
#         # else:
#         clean_df = clean_df.reset_index()
#         clean_df = clean_df.drop(columns=['index'])
#         #     clean_df = clean_df.rename(columns={'index':'idx'})
#         logging.info('L0 data pre-processed')

#     # Save processed data
#     if run_id is not None:
#         output_path = f'../data/L0/L0_{run_id}'
#         if output_format is not None:
#             if os.path.exists(output_path + output_format):
#                 logging.warning('Output file already exists. Overwriting...')
        
#         if output_format == '.feather':
#             clean_df.to_feather(output_path + '.feather')
#             logging.info('L0 .feather file created')
#         elif output_format == '.csv.zip':
#             clean_df.to_csv(output_path + '.csv.zip', index=False, compression='zip')
#             logging.info('L0 ".csv.zip" file created')
#         elif output_format == '.pkl.gzip':
#             clean_df.to_pickle(output_path + '.pkl.gzip', compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})
#             logging.info('L0 .pkl.gzip file created')
#         elif output_format is not None:
#             raise ValueError('Unsupported output file format')
#     elif output_format is not None:
#         raise ValueError('run_id is required for saving the output')
    

#     return clean_df


# def measure_execution_time(func, *args, **kwargs):
#     start_time = time.time()
#     result = func(*args, **kwargs)
#     end_time = time.time()
#     execution_time = end_time - start_time
#     return result, execution_time
