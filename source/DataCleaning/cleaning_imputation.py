import pandas as pd
import numpy as np
import re


PMCregex = r'PMC\d+'
PMIDregex = r'PMID:\s?\d+'

# The dictionary serves to convert publisher and journal names into the most
# common occurrences contained with the CSV file. For example converting
# abbreviated journal names to the full name. Note, these may not be accurate.
# The keys assume that text has been converted to lower case and then converted
# to title case.  In other words, the keys do not exist in the raw file until
# they have ran through the convert_journal_pub_names function.

PUB_JOURNAL_LOOKUP = {
    'Acs': 'American Chemical Society',
    'Acs (Amercian Chemical Society) Publications': 'American Chemical Society',
    'Asm': 'American Society for Microbiology',
    'Biomed Central Limited': 'BioMed Central',
    'Biomed Central Ltd':'BioMed Central',
    'Bmj':'British Medical Journal',
    'Bmj Group':'British Medical Journal',
    'Bmj Publishing Group':'British Medical Journal',
    'Bmj Publishing Group Ltd':'British Medical Journal',
    'Cadmus':'Cadmus Journal Services',
    'Cambridge Journals':'Cambridge University Press',
    'Cambridge Uni Press':'Cambridge University Press',
    'Cambridge Univ Press':'Cambridge University Press',
    'Cenveo Publisher Services/ASM JV1': 'Cenveo Publisher Services',
    'Cold Spring Habour Press':'Cold Spring Harbor Laboratory Press',
    'Cold Spring Harbor Press':'Cold Spring Harbor Laboratory Press',
    'Cold Spring Harbor Publications':'Cold Spring Harbor Laboratory Press',
    'Dartmouth Journals':'Dartmouth Journal Services',
    'Elseveier Science':'Elsevier',
    'Faseb':'Federation of American Societies for Experimental Biology',
    'Federation of American Societies for Experimental Biology (FASEB)':'Federation of American Societies for Experimental Biology',
    'Oup':'Oxford University Press',
    'Oxford Univ Pres':'Oxford University Press',
    'Oxford Journals':'Oxford University Press',
    'Oxford Journals (OUP)':'Oxford University Press',
    'Plos':'Public Library of Science',
    'Plos (Public Library of Science)':'Public Library of Science',
    'Royal Society': 'Royal Society of Chemistry',
    'Rsc':'Royal Society of Chemistry',
    'Rsc Publishing':'Royal Society of Chemistry',
    'Aids':'Journal of Acquired Immune Deficiency Syndromes',
    'Aids Journal':'Journal of Acquired Immune Deficiency Syndromes',
    'Aids UK':'Journal of Acquired Immune Deficiency Syndromes',
    'Jaids Journal of Acquired Immune Deficiency Syndromes':'Journal of Acquired Immune Deficiency Syndromes',
    'Journal of Acquired Immune Deficiency Syndroms (JAIDS)':'Journal of Acquired Immune Deficiency Syndromes',
    'Journal of AIDS':'Journal of Acquired Immune Deficiency Syndromes',
    'Plos 1':'Plos One',
    'Neuroimage':'NeuroImage'
}





def read_csv(file_path, encoding="ISO-8859-1"):
    """
    Simple function that invokes Pandas read_csv function.

    Args:
        file_path (str): Absolute file path to the CSV file to be read.

    Returns:
        data_frame: A Pandas data file with the contents of the CSV file.
    """
    data_frame = pd.read_csv(file_path, encoding)
    return data_frame


def retrieve_metadata(dataframe):
    """
     Extracts the shape, column names and number of null rows for
     each column of a Pandas data frame. This information is printed
     to standard ouput.

     Args:
        dataframe (dataframe): The Pandas data frame of interest

     Returns:
        Nothing
    """
    num_rows, num_columns = dataframe.shape
    columns = dataframe.columns.values.tolist()
    print('Data frame contains %d columns and %d rows' % (num_columns, num_rows))
    for elem in columns:
        print('column %s has %d empty rows'% (elem, sum(dataframe[elem].isnull())))


def _search_pmcid_regex(id_str):
    """
    Tries to match the string containing the journal ID against
    a PMCID regex. Returns the first instance of a match.

     Args:
        id_str (str): Cell contents of the ‘'PMID/PMCID' column.

     Returns:
       tuple: A tuple containing boolean and string containing
       the new ID.  The boolean may have a value of true indicating
       that a sub string matched the regular expression or false otherwise.
       If false the match str is empty.
    """
    result = False
    match_str = ''
    try:
        match_str = re.search(PMCregex, id_str).group()

        result = True
    except AttributeError:
        pass
    return result, match_str


def _search_pmid_regex(id_str):
    """
     Tries to match the string containing the journal ID against
    a Pub Med regex. Returns the first instance of a match.

     Args:
        id_str (str): Cell contents of the ‘'PMID/PMCID' column.

     Returns:
       tuple: A tuple containing boolean and string containing
       the new ID.  The boolean may have a value of true indicating
       that a sub string matched the regular expression or false otherwise.
       If false the match str is empty.
    """
    result = False
    match_str = ''
    try:
        match_str = re.search(PMIDregex, id_str).group()
        match_str = match_str.replace(':', '')
        match_str = match_str.replace(' ', '')
        result = True
    except AttributeError:
        pass
    return result, match_str


def pmcid_pmid_converter(id_str):
    """
     Specialized data massaging function for Wellcome trust data.
     This function should only be run on the first column containing
     the label ‘PMC ID/PMID’.
     The following assumptions and conversions are carried out;

        1. Those id numbers starting with 3 are PMCID ids
        2. Those id numbers starting with 2 are PMID ids
        3. Those ids without prefix will have either PMIC or PMID appended
        4. NA values will be kept
        5. Any other text in a row that does not adherer to the above two
           forms will be dropped and replaced with NA.
        6. Cells that have both PMIC and PMID identifiers will have
           the 2nd id dropped. Only one id is necessary.

     Args:
        id_str (str): String instance containing a cells entire text.

     Returns:
        str: string instance containing either NA, PMCID* or PMID*
    """
    data = str(id_str)
    data = data.strip('\n')
    has_pmcid_match = False
    has_pmid_match = False
    match_pmcid_str = ''
    match_pmid_str = ''
    result_str = ''

    has_pmcid_match, match_pmcid_str = _search_pmcid_regex(data)
    has_pmid_match, match_pmid_str = _search_pmid_regex(data)

    if has_pmcid_match or has_pmid_match:
        if has_pmcid_match:
            result_str = match_pmcid_str
        elif has_pmid_match:
            result_str = match_pmid_str
    elif data.isdigit():
        if data[0] == '2':
            result_str = 'PMID' + data
        elif data[0] == '3':
            result_str = 'PMCID' + data
    else:
        result_str = 'NA'

    return result_str


def remove_money_symbols(string):
    """
    Removes the pound symbol and trailing dollar sign symbols from
    the Cost data. Resulting data is converted to type float.

    Args:
        string (str): string instance housing cell contents.

    Returns:
        float: Cost of journal article represented as a float
    """
    string = string.replace('$', '')
    return float(string[1:])


def remove_cost_outliers(cost, dataf, threshold):
    """
    Sets any coat values that are equal to or greater than threshold
    to the 90th  quartile for the coat column.
    Assumes that a Seres object is passed in as the dataf argument
    having a single column.

    Args:
        coat (float):           cost value to be capped.
        dataf (Pandas, Series): Series instance that statistic will be retrieved from.
        threshold (decimal):    Outlier threshold

    Returns:
        float: Resulting cost value
    """
    result = cost
    if cost >= threshold:
        result = dataf.quantile(.9).values[0]
    return result


def convert_journal_pub_names(string):
    """
    Converts journal and Publication names into title case and transforms string
    into a generalized name. These generalized terms are housed in a dictionary
    and exchanged for the most least used term or mis spelled text. For example
    changing ‘Neuroimgae’ into ‘NeuroImage’

    Args:
        string (str): Publication or Journal column cell contents.

    Returns:
        str: The converted cell contents as a string instance.
    """
    string = str(string)
    string = string.lower()
    result = string.title()
    if result in PUB_JOURNAL_LOOKUP:
        result = PUB_JOURNAL_LOOKUP.get(result)
    return result
