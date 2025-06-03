### lab_calc_functions.py

import pandas as pd
import numpy as np
import re

prefix_dict = {"p":1e-12, "n":1e-9, "μ":1e-6, "u":1e-6, "":1, "m":1e-3, "c":1e-2, "k":1e3}
acc_bases = ["mol", "M", "g", "gram", "L", "liter", "m", "meter"]


def split_it(x):
    
    '''Uses regex to add a space between value and units if there isn't a space. Searches for a digit next to a letter 
    or a percentage sign &, if found, adds a space between them'''
    
    x = str(x)
    regex = r'([0-9]+)(\%|[A-Za-z]|μ)'
    return re.sub(regex, r'\g<1> \g<2>', x)


def get_sig_figs (value):
    '''takes value with units as string in format "value units" and returns number of significant figures'''
    
    value = split_it(value)
    
    if "." not in value:
        sig = value.strip("0")
    else:
        sig = value.replace(".","").lstrip("0")
    return len(sig)


def prefix2base(value_units):
    '''Takes value with units as string in format "value units" and returns the value in the SI base unit.'''
    value, units = value_units.split()
    if len(units) == 1:
        prefix = ""
        base = units
    else:
        prefix = units[0]
        base = units[1:]
    sig_figs = get_sig_figs(value)
    sci_not = np.format_float_scientific((prefix_dict[prefix]*float(value)), unique=False, precision=sig_figs-1)
    return sci_not+" "+base  


def fix_sig_figs(value):
    '''work in progress - takes a number and formats it in scientific notation with the current number of sig figs. '''
    sig_figs = get_sig_figs(str(value))
    sci_not = np.format_float_scientific((prefix_dict[prefix]*float(value)), unique=False, precision=sig_figs-1)
    return sci_not


def check_units(val_list):
    '''Takes a list of values and determines whether you need to convert any units. If yes, tells you which. Returns True if the units match, False if at least some of the units don't match'''
    parts = {}
    for i, val in enumerate(val_list):
        if len(val.split()) < 2:
            print("No units provided for", val)
            break
        unit = val.split()[1]
        if unit in acc_bases:
            prefix = ""
            base = unit
        else:
            prefix = str(unit[0])
            base = str(unit[1:])

        if base in parts:
            if parts[base]!=[prefix]:
                parts[base].append(prefix)
        else: 
            parts[base]=[prefix]

    fix_these = []

    for key in parts:
        if len(parts[key]) > 1:
            fix_these.append(key)
    if len(fix_these) > 0:
        print("You need to fix", fix_these)
        return False
    else:
        print("Units match.")
        return True

# if all units are compatible as is, compute without converting to base units
# else, convert all to base units

def convert_if_needed(vals):
    ''' Takes a list of values. If all units are compatible as is, returns values as is. Else, converts all to base units and returns.'''
    if not check_units(vals):
        conv_vals = []
        for i, val in enumerate(vals):
            conv_vals.append(prefix2base(val))
        print("I fixed them for you & converted everything to base units.")
        return conv_vals
    else:
        print("No conversions needed.")
        return vals


# take a value and determine what's the best metric prefix based on size
# choose prefix that avoids decimals

def best_prefix(value):
    '''Takes a numeric value and determine what's the best metric prefix based on size - chooses prefix that avoids decimals'''
    if value >= 1e3:
        return "k"
    elif 1 <= value <= 1e3:
        return ""
    elif 1e-3 <= value <= 1:
        return "m"
    elif 1e-6 <= value <= 1e-3:
        return "μ"
    elif 1e-9 <= value <= 1e6:
        return "n"
    elif value <= 1e9:
        return "p"


# convert between different prefixes
# default to base units if to_unit not provided

def metric_conversion(value_units, to_unit="base"):
    '''takes value with units as string in format "value units" and converts to a different metrix prefix. 
    Defaults to base units if to_unit not provided. If the user uses "best" as to_unit, the output will be converted
    to a unit that avoids decimals.
    
    -------------------------------------
    examples:
    
    defaults to converting to base unit
    
    >>>>>> metric_conversion("0.010 μM")
    '1.0e-08 M'
    
    alternatively, provide desired unit
    
    >>>>>> metric_conversion("0.010 μM", "mM")
    '1e-05 mM'
    '''
    
    value_base = prefix2base(value_units)
    value = value_base.split()[0]
    base_units = value_base.split()[1]
    
    #if to_unit is base (or not defined by user), convert to base
    if to_unit == "base" or to_unit in acc_bases:
        #print(str(round(float(value), sig_figs)) + " " + to_unit)
        return str(value)+ " " + base_units
    
    # if to_unit is best, convert to a prefix that avoids decimals
    elif to_unit == "best":
        prefix = best_prefix(float(value))
        return str(float(value)/prefix_dict[prefix])+ " " + prefix + base_units

    # otherwise, convert to the desired unit given
    else:
        prefix = str(to_unit[0])
        #print(str(round(float(value)/prefix_dict[prefix], sig_figs)) + " " + to_unit)
        return str(float(value)/prefix_dict[prefix])+ " " + to_unit

def calc_V1 (c1, v2, c2):
    '''Given initial concentration, final volume, and final concentration (as numeric), calculate and return initial volume.'''
    return (v2*c2)/c1

def calc_V2 (c1, v1, c2):
    '''Given initial concentration, inital volume, and final concentration (as numeric), calculate and return final volume.'''
    return (c1*v1)/c2

def calc_c1 (v1, c2, v2):
    '''Given initial volume, final concentration, and final volume (as numeric), calculate and return final concentration.'''
    return (c2*v2)/v1

def calc_c2 (c1, v1, v2):
    '''Given initial concentration, initial volume, and final volume (as numeric), calculate and return final concentration.'''
    return (c1*v1)/v2


# adapted function to default to returning values in unit that avoids decimals
# if the user specifies best=False, the output will be in base units, 
# unless the user specifies desired units with the des_units argument
# if des_units are given, best will be automatically converted to false

def dilution_calc(c1=False, v1=False, c2=False, v2=False, best=True, des_units="base", print_results=True):
    '''Takes 3-4 strings corresponding to inital and final concentrations and volumes
    in format "value units". If given 3 values, calculates and returns the remaining value.
    If given 4 values, checks whether the math works out. 
    
    Defaults to returning values in unit that avoids decimals. If the user specifies best=False, 
    the output will be in base units, unless the user specifies desired units with the des_units 
    argument. If des_units are given, best will be automatically converted to false.
    
    If prints_results==True (default), results will be printed. If prints_results==False, results will 
    be returned but not printed'''
   
    #first check if units are provided
    vals_in = [c1, v1, c2, v2]
    vals_in = [x for x in vals_in if x!=False]
    unit_check = [len(str(x).split()) > 1 for x in vals_in]
    #to do: determine precision from input
    precision = 1
    if any(unit_check):
        if print_results==True:
            print ("has units")
        #convert to base units for calculations -> later convert to best (or user-specified) units
        vals = []
        for val in vals_in:
            vals.append(prefix2base(val))
    else:
        if print_results==True:
            print ("no units")

    # if all values are provided, check if the math is correct
    if len(vals)==4:
        if c1*v1==c2*v2:
            return "Good job! You did my work for me!"
        else:
            return "Check your math..."

    # if not enough input, let them know if they didn't provide enough data to calculate a missing value
    elif len(vals) < 3:
        return "Not enough input given. A minimum of values are required."

    # if 3 values given, determine which value is missing, and calculate it
    else:
        # if the user doesn't want best prefix determined, default to reporting in base units
        # the user can alternatively specify the units they want used with the "des_units" argument

        # if user specifies a desired unit, default best to False even if they don't include best=False
        if des_units != "base":
            best = False

        if best == False:
            if des_units == "base":
                prefix = ""
            else:
                prefix = des_units[0]

        if not c1:
            if print_results==True:
                print("c1 =")

            v1 = float(vals[0].split()[0])
            c2 = float(vals[1].split()[0])
            v2 = float(vals[2].split()[0])

            value = calc_c1(v1, c2, v2)

            if best==False:
                value = value/prefix_dict[prefix]
                # if in base, report in scientific notation
                if prefix == "":
                    return np.format_float_scientific(value, unique=False, precision=precision)  + " " +  vals[1].split()[1]
                else: 
                    return str(value) + " " + prefix + vals[1].split()[1]

            else:
                prefix = best_prefix(value)
                value = value/prefix_dict[prefix]
                return str(value) + " " + prefix + vals[1].split()[1]

        elif not v1:
            if print_results==True:
                print("v1 =")

            c1 = float(vals[0].split()[0])
            c2 = float(vals[1].split()[0])
            v2 = float(vals[2].split()[0])

            value = calc_V1(c1, v2, c2)

            if best==False:
                value = value/prefix_dict[prefix]
                # if in base, report in scientific notation
                if prefix == "":
                    return np.format_float_scientific(value, unique=False, precision=precision)  + " " +  vals[2].split()[1]
                else:
                    return str(value) + " " + prefix + vals[2].split()[1]

            else:
                prefix = best_prefix(value)
                value = value/prefix_dict[prefix]
                return str(value) + " " + prefix + vals[2].split()[1]

        elif not c2:
            if print_results==True:
                print("c2 =")

            c1 = float(vals[0].split()[0])
            v1 = float(vals[1].split()[0])
            v2 = float(vals[2].split()[0])

            value = calc_c2(c1, v1, v2)

            if best==False:
                value = value/prefix_dict[prefix]
                # if in base, report in scientific notation
                if prefix == "":
                    return np.format_float_scientific(value, unique=False, precision=precision)  + " " +  vals[0].split()[1]
                else: 
                    return str(value) + " " + prefix + vals[0].split()[1]

            else:
                prefix = best_prefix(value)
                value = value/prefix_dict[prefix]
                return str(value) + " " + prefix + vals[0].split()[1]

        elif not v2:
            if print_results==True:
                print("v2 =")

            c1 = float(vals[0].split()[0])
            v1 = float(vals[1].split()[0])
            c2 = float(vals[2].split()[0])

            value = calc_V2(c1, v1, c2)

            if best==False:
                value = value/prefix_dict[prefix]
                # if in base, report in scientific notation
                if prefix == "":
                    return np.format_float_scientific(value, unique=False, precision=precision)  + " " +  vals[1].split()[1]
                else: 
                    return str(value) + " " + prefix + vals[1].split()[1]

            else:
                prefix = best_prefix(value)
                value = value/prefix_dict[prefix]
                return str(value) + " " + prefix + vals[1].split()[1]

# given mw & desired concentrations & volumes in form "value units", calculate g needed 

# conc is in M, so is mol/L
## mol needed = conc * vol
# the g needed to get that many mol...
## mw is g/mol
### so g needed is mol needed * mw
#### so, conc * vol * mw

def calc_mass_needed (mw, desired_conc, desired_vol):
    '''Given mw & desired concentrations & volumes in form "value units", calculates mass needed. Takes molecular weight as number (no units), desired concentration and desired volume with units as strings in format "value units." 
    
    uses this logic:
        conc is in M, so is mol/L
        mol needed = conc * vol
        the g needed to get that many mol...
        mw is g/mol
        so g needed is mol needed * mw
        so, conc * vol * mw'''
    
    conc_units = prefix2base(desired_conc)
    conc = float(conc_units.split()[0])
    vol_units = prefix2base(desired_vol)
    vol = float(vol_units.split()[0])
    
    g_needed = conc*vol*mw
    
    if g_needed >= 1:
        return str(g_needed) + " " + "g" 

    else:
        return str(g_needed*1e3) + " " + "mg"



# given mass in form "mass units", calculate mols, automatically determining
# what the prefix should be

# if the user specifies best=False, the output will be in base units, 
# unless the user specifies desired units with the des_units argument
# if des_units are given, best will be automatically converted to False

# mw is g/mol
# mol = mass/mw

def mass2mol(mw, mass, best=True, des_unit="base"):
    '''
    Given mass as string in form "mass units", calculates mols, automatically determining what the prefix should be. If the user specifies best=False, the output will be in base units, unless the user specifies desired units with the des_units argument. If des_units are given, best will be automatically converted to False.

    Uses this logic:
        mw is g/mol
        mol = mass/mw

    ---------------
    example:
    >>>> mass2mol(121.14, "60.57 g")
    '500.0 mmol'

    '''
    mass_units = prefix2base(mass)
    mass = float(mass_units.split()[0])
    
    mol = mass/mw
    
    # if user specifies a desired unit, default best to False even if they don't include best=False
    if des_unit != "base":
        best = False

    if best == False:
        if des_unit == "base":
            return str(mol) + " " + "mol"
        else:
            prefix = des_unit[0]
            return str(mol/prefix_dict[prefix]) + " " + prefix + "mol"

    else: 
        if mol >= 1:
            return str(mol) + " " + "mol"
        elif 1e-3 <= mol <= 1:
            return str(mol*1e3) + " " + "mmol"
        elif 1e-6 <= mol <= 1e-3:
            return str(mol*1e6) + " " + "μmol"
        elif 1e-9 <= mol <= 1e6:
            return str(mol*1e9) + " " + "nmol"
        elif 1e-12 <= mol <= 1e9:
            return str(mol*1e12) + " " + "pmol"


# determine molarity of a solution given the mw, mass, & volume
# mw should be given as numeric, no units
# mass & volume should be given as strings, with space between value and units (i.e., "value units")

# we want mol/L
# we know mw (g/mol), mass (which we will convert to g), and volume (which we will convert to L)
## so, basically we know g/L
### to get to mol/L, we need to convert g to mol, which we can do by dividing by the mw
#### if we do (g/L)/(g/mol), g cancel out, leaving us with mol/L

# the most challenging part is just dealing with formatting & units
# we will do the calculations with values in the base units (M, g, L)
# we can then convert to a more practical unit based on how big the values are
## alternatively, user can specify desired units

def weightvol2molar (mw, mass, vol, des_units="best"):
    '''Determines molarity of a solution given the mw, mass, & volume. mw should be given as numeric, no units. mass & volume should be given as strings, with space between value and units (i.e., "value units"). Molarity  will be given with a prefix that avoids decimals unless user specifies units with "des_units"
    
    --------------
    example:
    >>>> weightvol2molar(121.14, "60.57 g", "500 mL", "mM")
    '1000.0mM'
    '''
    
    mass_units = prefix2base(mass)
    mass = float(mass_units.split()[0])
    vol_units = prefix2base(vol)
    vol = float(vol_units.split()[0])
    
    molarity = (mass/vol)/(mw)
    
    if des_units == "best":
        # convert to a prefix that avoids decimals
        prefix = best_prefix(molarity)
        molarity = molarity/prefix_dict[prefix]
        return str(molarity) + " " + prefix + "M"
    
    else:
        prefix = des_units[0]
        molarity = molarity/prefix_dict[prefix]
        return str(molarity) + " " + prefix + "M"
    
# determine mass/mL of a solution given the mw & molarity
# mw should be given as numeric, no units
# molarity should be given as a string, with space between value and units (i.e., "value units")

# we want mass/vol (g/L in base units)
# we know mw (g/mol) and molarity (mol/L)
### to get to g/L, we just need to multiply these, and mol will cancel out
# we can then convert to a more practical unit based on how big the values are
## alternatively, user can specify desired units

def molar2weightvol (mw, molarity, des_units="best"):
    '''Determines weight/volume concentration of a solution given the mw & molarity. mw should be given as numeric, no units. molarity should be given as a string, with space between value and units (i.e., "value units"). Weight volume  will be given with a prefix that avoids decimals unless user specifies units with "des_units"
    '''
    
    molarity_units = prefix2base(molarity)
    molarity = float(molarity_units.split()[0])
    
    # first get mass/vol in g/L, then divide by 1000 to get g/mL
    mass_vol = (mw*molarity)/1000
    
    # then convert the g to more convenient units
    
    if des_units == "best":
        # convert to a prefix that avoids decimals
        prefix = best_prefix(mass_vol)
        mass_vol = mass_vol/prefix_dict[prefix]
        return str(mass_vol) + " " + prefix + "g/mL"
    
    else:
        prefix = des_units[0]
        mass_vol = mass_vol/prefix_dict[prefix]
        return str(mass_vol) + " " + prefix + "g/mL"


# given mw & desired concentration & mass in form "value units", calculate volume to dissolve in

# conc is in M, and M = mol/L
# we can rearrange this to L = mol/M

# we have M, so we just need to calculate mol
## this will be done by calling the mass2mol function, which does the following calculation
## we can calculate mol from mw & mass
## mw is in g/moL 

## mw = g/moL
## we can rearrange this to
## mol = g/mw

# then we just have to divide the mol by the concentration to get L
# L = (mol/M)

# and then we can format it to return in a convenient form

def calc_vol_needed (mw, mass, desired_conc):
    '''Given mw & desired concentration & mass as strings in form "value units", 
    calculates the volume to dissolve in and returns the volume in mL if < 1L and 
    L if > 1L.'''
    mass_units = prefix2base(mass)
    mass = float(mass_units.split()[0])

    conc_units = prefix2base(desired_conc)
    conc = float(conc_units.split()[0])

    mol = float(mass2mol(mw, mass_units, best=False).split()[0])

    vol_needed = mol/conc

    if vol_needed >= 1:
        return str(vol_needed) + " " + "L"

    else:
        return str(vol_needed*1e3) + " " + "mL"
    
def dilute_stock(compound, des_conc, des_vol, mol_weights, stock_conc):
    '''takes information about a stock compound from mol_weights and stock_conc dictionaries (with
    "compound" as the corresponding key), and desired concentration (in mg/mL) & volume & returns
    the volume of stock and volume of solvent needed'''
    
    c1_stock = weightvol2molar(mol_weights[compound], stock_conc[compound], "1 mL")
    v1_stock = dilution_calc(c1=c1_stock, c2=des_conc, v2=des_vol)
    v_sol = float(des_vol.split()[0]) - float(v1_stock.split()[0])
    return v1_stock, v_sol
    
# determine doubling time (e.g. for cell culture)

def doubling_time(c1, c2, t):
    '''takes unitless values for initial concentration, final concentration, and time and returns doubling time. all numeric'''
    dt = (t*log(2))/(log(c2)-log(c1)) 
    return round(dt, 1)


def buffer_brewer(reagent_info_dict, v_desired, x_factor=1):
    
    '''Calculates volumes of stock solutions needed for desired final solution. 
    Returns a dataframe with initial concentrations, final concentrations, & needed
    input volumes. Also prints out a formatted list of what to prepare. Input can
    be given as a dictionary, dataframe, or csv file.
    
    Dictionary should be in format: {"reagent 1": ["initial concentration", "final concentration"], 
    "reagent 2": ["initial concentration", "final concentration"], ...}
    
    e.g. {"Tris pH 8.0": ["1 M", "50 mM"], "NaCl": ["5 M", "150 mM"], 
                      "Triton X-100": ["100 %", "1 %"], "sodium deoxycholate": ["10 %", "0.5 %"], 
                      "SDS": ["10 %", "0.1 %"]}
    
    DataFrame or csv should have reagent names as index, then initial concentration as first column
    and final concentration as second column.
    
    x_factor can be specified in order to make stock solutions. It defaults to 1 if not given
    
    # work in progress: If there isn't a space between the value and the unit, add one
    
    columns are renamed in case the input names are different
    
    '''
    
    # create dataframe if needed
    if isinstance(reagent_info_dict, pd.DataFrame):
        reagent_info = reagent_info_dict
        reagent_info.columns = ["C_initial", "C_final"]
    
    elif isinstance(reagent_info_dict, dict):
        # collect data into pandas dataframe
        reagent_info = pd.DataFrame(reagent_info_dict).T
        reagent_info.columns = ["C_initial", "C_final"]
    
    else:
        reagent_info = pd.read_csv(reagent_info_dict, index_col=0)
        reagent_info.columns = ["C_initial", "C_final"]
    
    # if they give their info as a csv
    
    ## TO DO
    #### use regex to add a space between value and units if there isn't a space
    #### search for a digit next to a letter or a percentage sign
    #### if found, add a space between them
    
    # parse v_desired for numeric value and units
    v_final = float(v_desired.split()[0])
    v_units = v_desired.split()[1]
    
    # adjust desired concentrations for x-factor given (defaults to 1)
    
    reagent_info["C_final"] = reagent_info["C_final"].apply(lambda x: str(float(x.split()[0])*x_factor) + " " + x.split()[1])
    
    # perform calculations to find initial volumes needed
    reagent_info["V_initial"] = reagent_info.apply(lambda row: dilution_calc(
        c1= row["C_initial"], c2= row["C_final"], v1=False, v2= v_desired, print_results=False), axis=1)
    
    # determine volume of water needed to dilute to final volume
    
    ## find combined volume of other components - convert the values to floats, convert them
    ## to base units (L), and sum them together
    non_water_vol = sum(reagent_info["V_initial"].apply(prefix2base
                                 ).apply(lambda x: x.split()[0]
                                        ).astype(float))
    
    ## subtract the non-water volume from the total volume to find the water volume
    v_final_base = prefix2base(v_desired)
    water_vol = str(float(v_final_base.split()[0]) - non_water_vol) + " L"
    water_vol = metric_conversion(water_vol, to_unit="best")
    
    # add water to dataframe
    reagent_info.loc["water"] = [np.nan, np.nan, water_vol]
    
    # format decimals displayed
    reagent_info["V_initial"] = reagent_info["V_initial"].apply(lambda x: '{0:.2f}'.format(float(x.split()[0])) 
                                                      + " " + x.split()[1])
    
    # print it nicely formatted
    formatted_starts = reagent_info.apply(lambda row: row.name + 
                                                          " (" + str(row["C_initial"]) + ")", axis=1)
    to_prepare = formatted_starts + ": " + reagent_info["V_initial"].astype(str)
   
    if x_factor != 1:
        print("To prepare " + v_desired + " " + str(x_factor) + " X:")
    else:
        print("To prepare " + v_desired + ":")
    print(to_prepare.to_string(index=False))
    
    return reagent_info


def buffer_from_solids(mol_weights, recipe_dict, desired_vol):
    '''Determines how much of each component to dissolve to get a desired solution. 
    
        Takes as imput:
            mol_weights: dictionary of molecular weights
            
            recipe_dict: a dictionary of the buffer recipe containing reagent name as the key and the 
            desired concentration as the value (in the form "value units")
            
            desired_vol: desired volume in the form "value units" '''

    
    for component in recipe_dict:
        print(component, calc_mass_needed(mol_weights[component], recipe_dict[component], desired_vol))