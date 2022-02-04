'''
Hydrus input file generator for running ssat v. saturate paste extract
Created by A.J. Brown and Cullen McGovern on 10/11/18 in python 3.6.4

To use:
place this module into folder where 'H1D_UNSC.exe' is located, and place
input file templates into the project folder where inputs and outputs are
stored when hydrus is ran (i.e. the 'Projects' folder).

Update Log:
12/5/18
    -Merged all functions into class.
    -Hyrdrus runs over all iterations, but ECe is not changing (as it should).
     Not sure why.
12/12/18
    -Hydrus completes one iteration successfully, but cannot write output files
    -HydrusRun_ver2 was created to attempt running from H1D_UNSC.EXE directory
    -runs one iteration successfully, and writes complete output files, but
     cannot run the next iteration.
12/18/18
    -Hydrus completes all iterations successfully! Now we need to check the
     quality of the data by plotting some things:
    -ECe v. Gypsum content with colorization for water contents for one soil
    -ECe v. Gypsum content with colorization for water contents for all soils
    New Efficiency Issues:
    -Running graph_output() results in having to run final_data() 4 separate
     times because of the way it is called with self.final_data() in the
     function. There is probably a way to prevent this from happening by saving
     the dataframe after running final_data() just once.
12/20/18
    -Moving a little backwards. We decided it wasn't good to just select the
     first EC from the Equil.out file and instead import all of it, in case in
     the future someone might want to look at the entire soil profile at
     different times. This resulted in changing the _run_batch() function, but
     it's not working properly.
    -_run_batch() creates a giant DataFrame, but should only produce one of 1176
     (294 rows from input * 4 rows from Equil.out = 1,176 rows). Instead we have
     a massive dataframe with 43,932 rows, but we're unsure why.
01/02/2019
    -fixed merge issue. More than one index was needed to merge correctly.
     So instead of merging on just 'soil_type', gypsum and VWC were also merged
     on.
    -consolidated ssat and paste ion values into single columns to allow hydrus
     to run each sample separately.  This will make hydrus simulate a paste, and
     then simulate the field condition without mising ions from both.
    -implemented graph_output() function with customizable x and y values
    -need to improve graph_output() to include separate line fits for paste and
     ssat points
02/11/19
    -need to figure out why NaN's are produced in input dataframe
        -still unknown
    -need to optimize gypsum so that the 1:1 fit is better
        -this might not be the case, I might just take observed values from the
         CSU soils lab for both CaSO4 and CaCO3
    -then I need to compare the same sample at different water contents
        -completed, added this functionality in self._inputs()
02/12/19
    -need to update self.graph_output() so that we can graph different moistures
        -think that faceting would be good?
        -can't figure out why slicing my group isn't working in the graph_output
02/18/19
    -fixed the graphing a bit
        -lost colorization
        -still can't condition by shape
        -can't separate moisture levels
    -maybe I don't need to group anything? I can probably just graph the df
    -DONE: try to include CO2 into datafram as input for paste vs. ssat
        -ambient CO2 (410ppm), ssat CO2 (410*1.05 = 430.5ppm)
        -CO2 is lower in lab water than in soil water
        -Im confused what hydrus CO2 is: ambient or solution?
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

WORKING_DIR = (
    'C:/Users/Public/Documents/PC-Progress/Hydrus-1D 4.xx/Projects/Gypsum')
HYDRUS_DIR = '.'

class HydrusRun:
    def __init__(self, projectDir=WORKING_DIR, hydrusDir=HYDRUS_DIR,
        soilFile='Gapon Output 2019-02-15.csv',
        selectorTemplate='SELECTOR0.IN',
        profileTemplate='PROFILE_2node_CO2.DAT', #changed for only 2 soil nodes
        executable='H1D_UNSC.EXE',
        outputFile='Equil.out',
        gypsumLimit=0, timeout=60):
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
        _inputs() row# =
        (7 soils*2 soil types(paste & ssat)*21 gypsum levels*2 moisutres)
        = 588 rows
        '''
        df = pd.concat([
            self._build_gypsum(r) for _, r in self.soilData.iterrows()])
        cols = df.columns
        #data into long format and gets rid of 'ssat' and 'paste' prefixes
        paste = (
            df[[c for c in cols if not c.lower().startswith('ssat')]]
            .rename(columns={c: c[5:] for c in cols if c.lower()
            .startswith('paste')}).assign(paste=True, CO2=0.000431))
        ssat = (
            df[[c for c in cols if not c.lower().startswith('paste')]]
            .rename(columns={c: c[4:] for c in cols if c.lower()
            .startswith('ssat')}).assign(paste=False, CO2=0.00041))
        #Adds the 2 levels of moisture content
        p = paste.assign(field='lab').append(paste.assign(VWC=ssat.VWC, field='field'))
        s = ssat.assign(field='field').append(ssat.assign(VWC=paste.VWC, field='lab'))
        #return p,s # For Debugging
        return pd.concat([p, s]).reset_index(drop=True)

    def _run_batch(self):
        '''
        Function to create input files SELECTOR.IN and PROFLE.DAT, and then run
        hydrus for each row from the _inputs dataframe. The output is then
        collected from the equil.out file and combined into a master dataframe.
        for example 28 input rows * 4 output rows per input = 112 master rows
        '''
        inputs = []
        outputs = []
        for idx, data in self._inputs().iterrows():
            inputs.append(data)
            run = data.to_dict()
            self._build_selector(run)
            self._build_profile(run)
            row = self._run().assign(
            soil_type=data.soil_type, gypsum=data.gypsum, VWC=data.VWC,
            paste=data.paste, field=data.field)
            outputs.append(row)
            #print('Iteration No: ', idx + 1, os.linesep)
        inputs = pd.DataFrame(inputs)
        outputs = (
            pd.concat(outputs)
            .rename(columns={'i':'hydrus_node', 'x':'hydrus_depth',
                'SAR':'hydrus_SAR'}))
        merged = (
            outputs.merge(
                inputs, on=['soil_type', 'gypsum', 'VWC', 'paste', 'field'],
                how='left'))
            #.set_index('soil_type'))
        #return inputs, outputs, merged  # For debugging
        return merged.reset_index(drop=False)

    def graph_output(self, sat_paste='ECe', ssat='ECpw', y='EC',
                        legend=True, line=2, time=0, depth=0):
        #call whole dataframe
        df = self._run_batch()
        #slices data so that the soil is at desired time and depth in hydrus
        df = df.loc[(df.time == time) & (df.hydrus_depth == depth)]
        df_paste = df[(df.paste==True) & (df.field=='lab')]
        df_ssat = df[df.paste==False  & (df.field=='field')]
        #inititalizes empty plot
        fig, ax = plt.subplots()
        ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling

        #Colorizes by soil type (different color is given to each soil)
        n = len(df.index.unique())
        cmap = plt.cm.get_cmap('hsv', n+1)

        #plots scatter points
        Paste = plt.scatter(df_paste[sat_paste], df_paste[y],
                    c='orange', marker='X')
        SSAT = plt.scatter(df_ssat[ssat], df_ssat[y], c='cyan',
                    marker='o')



        #determines if line of best fit (LOBF) is shown
        if line == 0:
            plt.plot()

        elif line == 1: #this is to create a LOBF for all points
            x = pd.concat([df[sat_paste], df[ssat]])
            y1 = pd.concat([df[y], df[y]])
            b, m = polyfit(list(x), list(y1), 1)
            plt.plot(x, b + m * x, '-')

        elif line == 2: #this is to create 2 lines for ssat AND sat_paste
            x, x1 = df_ssat[ssat], df_paste[sat_paste]
            y1, y2 = df_ssat[y], df_paste[y]
            b, m = polyfit(list(x), list(y1), 1)
            b1, m1 = polyfit(list(x1), list(y2), 1)
            plt.plot(x, b + m * x, '-', color='black')
            plt.plot(x1, b1 + m1 * x1, '-.', color='black')

        else: #this option DOES NOT graph a LOBF
            raise Exception('Please choose a valid line count: 0, 1, or 2')

        #determines if legend is shown
        if legend == True:
            ax.legend([Paste, SSAT],
                      [sat_paste, ssat, 'SSAT LOBF', 'Paste LOBF'])

        #plot labelling
        plt.axis('equal')
        plt.suptitle('Graph of ' + sat_paste + ' and ' + ssat + ' v. ' + y)
        if sat_paste != ssat:
            plt.xlabel(sat_paste + ' and ' + ssat)
        else:
            plt.xlabel(ssat)
        plt.ylabel(y)
        plt.show()
        '''
        It is worth noting here that SmeqL has a HUGE correlation with EC.
        Iterating over SmeqL instead of gypsum gives a VERY STRONG linear
        correlation. This could suggest a possible correction ion for ECe to
        ECpw (and thus EC of a saturated field).

        Plotting EC's from hydrus at different VWC seems to suggests a similar
        non-linear shift.  Try plotting lab ECe with lab SSAT and see. Then try
        field ECe with field SSAT.
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
            # join lines, encode UTF-8, and read it all into file-like object
            # for read_fwf
            df = pd.read_fwf(
                BytesIO(linesep.join(lines[start:stop]).encode()),
                skiprows=[1]) # skip the units row
            # add time column (why was this so hard for the authors?
            # Obviously a totally made up format is better...)
            df['time'] = time
            # add to results
            res.append(df)
        return pd.concat(res)

    def _build_gypsum(self, soil):
        #creates dataframe with a particular soil and increments of gypsum conc.
        res = []
        for n in range(self.gypLim):
            tmp = soil.copy()
            tmp['gypsum'] = n
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
