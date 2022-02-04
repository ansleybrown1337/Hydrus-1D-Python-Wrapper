'''
Hydrus input file generator for running ssat v. saturate paste extract
Created by A.J. Brown and Cullen McGovern on 10/11/18 in python 3.6.4

To use: place this module into folder where 'H1D_UNSC.exe' is located and place
input file templates into the project forlder where inputs and outputs are
stored when hydrus is ran.


Update Log:
10/11/18
Created
12/5/18
Merged all functions into class.
Hyrdrus runs over all iterations, but ECe is not changing (as it should).
Not sure why.
12/12/18
Hydrus completes one iteration successfully, but cannot write output files
HydrusRun_ver2 was created to attempt running from H1D_UNSC.EXE directory
-runs one iteration successfully, and writes complete output files, but cannot
 run the next iteration.
12/18/18
Hydrus completes all iterations successfully! Now we need to check the quality
of the data by plotting some things:
-ECe v. Gypsum content with colorization for water contents for one soil
-ECe v. Gypsum content with colorization for water contents for all soils
New Efficiency Issues:
-running graph_output() results in having to run final_data() 4 separate times
 because of the way it is called with self.final_data() in the function. There
 is probably a way to prevent this from happening by saving the dataframe
 after running final_data() just once.
12/20/18
-Moving a little backwards. We decided it wasn't good to just select the first
 EC from the Equil.out file and instead import all of it, in case in the future
 someone might want to look at the entire soil profile at different times. This
 resulted in changing the _run_batch() function, but it's not working properly
-_run_batch() creates a giant DataFrame, but should only produce one of 1176
 (294 rows from input * 4 rows from Equil.out = 1,176 rows). Instead we have a
 massive dataframe with 43,932 rows, but we're unsure why.
01/02/2019
-consolidated ssat and paste ion values into single columns to allow hydrus
 to run each sample separately.  This will make hydrus simulate a paste, and
 then simulate the field condition without mising ions from both.
-implemented graph_output() function with customizable x and y values
'''
#libraries
import os
import sys
import time
import datetime
import subprocess

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
from pathlib import Path
from io import BytesIO
from os import linesep

WORKING_DIR = 'C:/Users/Public/Documents/PC-Progress/Hydrus-1D 4.xx/Projects/Gypsum'
HYDRUS_DIR = '.'

class HydrusRun:
    def __init__(self, projectDir=WORKING_DIR, hydrusDir=HYDRUS_DIR,
        soilFile='Gapon Output 2019-01-03.csv',
        selectorTemplate='SELECTOR0 - SmeqL.IN',
        profileTemplate='PROFILE_2node.DAT', #changed for only 2 soil nodes
        executable='H1D_UNSC.EXE',
        outputFile='Equil.out',
        gypsumLimit=50, timeout=60):
        '''
        This is super cool docstringage
        '''
        self.timeout = timeout
        self.executable = executable
        self.gypLim = gypsumLimit + 1
        self.hydrusDir = Path(hydrusDir).absolute()
        self.projectDir = projectDir = Path(projectDir)
        self.selectorPath = selectorPath = projectDir.joinpath(selectorTemplate)
        self.profilePath = profilePath = projectDir.joinpath(profileTemplate)
        self.soilData = pd.read_csv(projectDir.joinpath(soilFile))
        self.outputFile = outputFile
        with open(selectorPath) as f:
            self.selectorText = f.read()
        with open(profilePath) as f:
            self.profileText = f.read()
        #file neded for hydrus to run outside of the GUI.
        #Does not exist by default.
        #Should exist in the same directory as the python module
        with open(self.hydrusDir.joinpath('LEVEL_01.DIR'), 'w') as f:
            f.write(str(self.projectDir))


    def _inputs(self):
        '''
        Function that takes input "Gapon Output 2018-12-05.csv" and for each
        row (i.e. soil type), generates 21 increments of gypsum from 0-20 meq/kg
        at 2 water contents (i.e. field VWC and saturated paste VWC). This
        totals 21 rows added for each original row of the input. For example:
        if original row# = 7, then:
        _inputs() row# = (7 soils * 2 water contents * 21 gypsum levels) = 294
        '''
        df = pd.concat([
            self._build_Sulfur_meqL(r) for _, r in self.soilData.iterrows()])
        cols = df.columns
        paste = (
            df[[c for c in cols if not c.lower().startswith('ssat')]]
            .rename(columns={c: c[5:] for c in cols if c.lower().startswith('paste')})
            .assign(paste=True))
        ssat = (
            df[[c for c in cols if not c.lower().startswith('paste')]]
            .rename(columns={c: c[4:] for c in cols if c.lower().startswith('ssat')})
            .assign(paste=False))
        return pd.concat([paste, ssat]).reset_index()

    def _run_batch(self):
        '''
        Function to create input files SELECTOR.IN and PROFLE.DAT, and then run
        hydrus for each row from the _inputs dataframe.
        '''
        inputs = []
        outputs = []
        for idx, data in self._inputs().iterrows():
            inputs.append(data)
            run = data.to_dict()
            self._build_selector(run)
            self._build_profile(run)
            row = self._run().assign(
            soil_type=data.soil_type, Sulfur_meqL=data.Sulfur_meqL, VWC=data.VWC)
            outputs.append(row)
            #print('Iteration No: ', idx + 1, os.linesep)
        inputs = pd.DataFrame(inputs)
        outputs = (
            pd.concat(outputs)
            .rename(columns={'i':'hydrus_node', 'x':'hydrus_depth',
                'SAR':'hydrus_SAR'}))
        merged = (
            outputs.merge(
                inputs, on=['soil_type', 'Sulfur_meqL', 'VWC'], how='left')
            .set_index('soil_type'))
        return merged

    def graph_output(self, sat_paste='Sulfur_meqL', ssat='Sulfur_meqL', y='EC',
                        legend = True):
        df = self._run_batch()
        groups = (
            df.loc[(df.time == 0) & (df.hydrus_depth == 0)]
            .groupby(['soil_type', 'paste']))
        fig, ax = plt.subplots()
        ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
        x = pd.concat([df[sat_paste], df[ssat]])
        y1 = pd.concat([df[y], df[y]])
        b, m = polyfit(list(x), list(y1), 1)
        plt.plot(x, b + m * x, '-')
        for keys, group in groups:
            soil_type, paste = keys
            if paste == True:
                ax.scatter(group[sat_paste], group[y], c='blue', label=keys,
                            marker='.')
            else:
                ax.scatter(group[ssat], group[y], c='red', label=keys,
                            marker='+')
        if legend == True:
            ax.legend()
        plt.suptitle('Graph of ' + sat_paste + ' and ' + ssat + ' v. ' + y)
        plt.xlabel(sat_paste + ' and ' + ssat)
        plt.ylabel(y)
        plt.show()
        '''
        It is worth noting here that SmeqL has a HUGE correlation with EC.
        Iterating over SmeqL instead of gypsum gives a VERY STRONG linear
        correlation. This could suggest a possible correction ion for ECe to
        ECpw (and thus EC of a saturated field).
        '''

    def _run(self):
        '''
        Find's hydrus from the installation directory, executes a model run, and
        collects saturated paste electrical conductivity (ECe) from the output
        file, Equil.out, at the top of the soil profile (in this case, ECe is
        the same throughout the entire profile so the first row will suffice).
        '''
        cmd = str(self.hydrusDir.joinpath(self.executable))
        run = subprocess.run(cmd, check=True,# timeout=self.timeout,
            stdout=subprocess.PIPE)
        return self._read_equil()

    def _read_equil(self):
        with open(os.path.join(self.projectDir,self.outputFile),'r') as f:
            lines = f.readlines()
        # drop all the newline charachters
        lines = [l.rstrip() for l in lines]
        # get rid of empty lines and ridiculous 'end' line
        lines = [l for l in lines if l and l != 'end']
        # find the lines that start with 'Time'
        times = [l for l in lines if l.startswith('Time')]
        # collect indices of times
        idxs = [lines.index(t) for t in times]
        # adjust indices to account for dropped times
        idxs = [i - 1 if i > 0 else 0 for i in idxs]
        # drop times
        lines = [l for l in lines if not l in times]
        # Now convert times to actual numbers
        times = [float(t.split(' ')[-1]) for t in times]
        # collect frames
        res = []
        for i, time in enumerate(times):
            # start where the time was
            start = idxs[i]
            # stop at the next time
            try:
                stop = idxs[i+1]
            # if we've gone past the end of idxs
            except IndexError:
                # list[i:None] is equivalent to list[i:]
                stop = None
            # join lines, encode UTF-8, and read it all into file-like object for read_fwf
            df = pd.read_fwf(
                BytesIO(linesep.join(lines[start:stop]).encode()),
                skiprows=[1]) # skip the units row
            # add time column (why was this so hard for the authors? Obviously a totally made up format is better...)
            df['time'] = time
            # add to results
            res.append(df)
        return pd.concat(res)

    def _build_Sulfur_meqL(self, soil):
        #creates dataframe with a particular soil and increments of Sulfur_meqL conc.
        res = []
        for n in range(self.gypLim):
            tmp = soil.copy()
            tmp['Sulfur_meqL'] = n
            res.append(tmp)
        return pd.DataFrame(res)

    def _build_selector(self, run):
        #bulds "SELECTOR.IN" file from template named under __init__()
        with open(self.projectDir.joinpath('SELECTOR.IN'), 'w') as f:
            f.write(self.selectorText.format(**run))

    def _build_profile(self, run):
        #bulds "PROFILE.DAT" file from template named under __init__()
        with open(self.projectDir.joinpath('PROFILE.DAT'), 'w') as f:
            f.write(self.profileText.format(**run))
