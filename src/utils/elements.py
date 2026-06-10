import re, os
import pandas as pd

_DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')

def sort_analytes(method, analytes, order = 'd'):
    """Sort the analyte list

    Sorting the analyte list can make data selection easier, or improve the pattern of correlations and PCA vectors.

    Parameters
    ----------
    method : str
        Method used for sorting.  Options include ``'alphabetical'``, ``'atomic number'``, ``'mass'``, ``'compatibility'``, and ``'radius'``.
    analytes : list
        List of analytes to sort
    order : str, optional
        Sets order as ascending (``'a'``) or decending (``'d'``), by default 'd'

    Returns
    -------
    list
        Sorted analyte list.
    """
    # Extract element symbols and any mass numbers if present
    parsed_analytes = []
    for analyte in analytes:
        # Extracts the element symbol and mass if available (e.g., "Al27" -> ("Al", 27))
        match = re.match(r"([A-Za-z]+)(\d*)", analyte)
        element_symbol = match.group(1) if match else analyte
        mass_number = int(match.group(2)) if match.group(2) else None
        parsed_analytes.append((element_symbol, mass_number))
    
    # Convert to DataFrame for easier manipulation
    df_analytes = pd.DataFrame(parsed_analytes, columns=['element_symbol', 'mass'])

    sort_data = pd.read_excel(os.path.join(_DATA_DIR, 'element_info.xlsx'))

    # Merge with sort_data for additional information
    df_analytes = df_analytes.merge(sort_data, on='element_symbol', how='left')

    # Sort based on the selected method
    match method:
        case 'alphabetical':
            df_analytes.sort_values(by='element_symbol', ascending=True, inplace=True)
        case 'atomic number':
            df_analytes.sort_values(by='atomic_number', ascending=True, inplace=True)
        case 'mass':
            # Use provided mass or average mass if not available
            df_analytes['computed_mass'] = df_analytes['mass'].fillna(df_analytes['average_mass'])
            df_analytes.sort_values(by='computed_mass', ascending=True, inplace=True)
        case 'compatibility':
            df_analytes.sort_values(by='order', ascending=False, inplace=True)
        case 'radius':
            df_analytes.sort_values(by='radius1', ascending=True, inplace=True)
        
    analytes = df_analytes['element_symbol'] + df_analytes['mass'].astype(str)
    # Return the sorted list of analytes as (symbol, mass) tuples
    return analytes.to_list()