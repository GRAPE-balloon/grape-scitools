contains software tools to process GRAPE data collected with PETsys system. 

As of 2024/06/11 the data folders are here for the purposes of running the code with required directory structure - data is stored elsewhere. 

# GRAPE Data Management 

## Processing Levels

**The Gamma Ray Polarimeter Experiment (GRAPE)** data products are processed at various levels ranging from Level 0 to Level 5. Raw data refers to unprocessed data acquired from the daq system. Level 0 products are raw data at full instrument resolution, with any and all acquisition artifacts removed. At higher levels, the data are processed and converted into more useful parameters and formats. As of 12/14/2023 data have not been properly formatted to meet FITs standards.

### Data Level Description

Level | Description
---------|----------
 L0 | Time stamped (flight data), reconstructed (if from binary), unprocessed instrument data at full resolution, with any and all acquisition artifacts removed (e.g., duplicate data, channel glitches, electronic anomalies), and raw 'energy' values are represented in 'adc' values.  
 L1 | L0 data annotated with ancillary information, including variables that provide identifying information about each detector. This can include but is not limited to scintillator type, location within the instrument, location within the column, associated temperatures. 
 L2 | Processed L1 data. 
 L2_1-N | L1 N-hit data
 <!-- L3A-Z | L2A-Z data calibrated to energy units (keV)
 L4A-Z | L3A-Z processed data discriminated/thresholds applied
 L5  | output or results from analyses of lower-level data -->

### Ancillary Data Description

 Level | Description
---------|----------
 Ancillary Data 1 | Data file containing variables representing identifying information about each detector such as scintillator type, location within the instrument, location within the column, asic, pins, gain identifiers. calibration coefficients computed and appended but not applied to L0 data.
 Ancillary Data 2 | Data file containing variables representing identifying information about each detector such as scintillator type, location within the instrument, and location within the column and computed calibration coefficients.
 Ancillary Data 3 | Data file produced with Petsys "convert raw to raw" method to extract debugging parameters: charge integration information.

###

Term | Description
---------|----------
Recorded data directory | A file directory containing all the recorded run data consisting a collection of files and metadata used by the Petsys extraction routine to compile 'raw' data from each data run.    
Raw data | This is the raw data we see - for each run this corresponds to a single data file in either binary, text, or root format. The raw data files contain information corresponding to time measurements, charge measurements, and detector ID.
Edited data | Data that has been corrected for errors and indexed. Flight data are tagged with a unix epoch (UTC). 
