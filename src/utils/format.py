import re
import numpy as np

def oround(val, order=2, toward=None):
    """Rounds a single of value to n digits

    Round a single number to the first n digits.

    Parameters
    ----------
    val : float
        number to round
    order : int, optional
        number of orders to retain, by default 2
    toward : int, optional
        direction for rounding, ``0`` for down, ``1`` for up, and ``None`` for nearest, by default None

    Returns
    -------
    float
        rounded matrix values
    """
    if val == 0:
        return 0

    power = np.floor(np.log10(abs(val)))
    if toward == 0:
        return np.floor(val / 10**(power-order)) * 10**(power - order)
    elif toward == 1:
        return np.ceil(val / 10**(power-order)) * 10**(power - order)
    elif toward is None:
        return np.round(val / 10**(power-order)) * 10**(power - order)
    else:
        raise ValueError("Valid values of toward may be 0, 1, or None.")

def oround_matrix(val, order=2, toward=None):
    """Rounds a matrix of values to n digits

    Rounds all numbers in a matrix to the first n digits.

    Parameters
    ----------
    val : float
        number to round
    order : int, optional
        number of orders to retain, by default 2
    toward : int, optional
        direction for rounding, ``0`` for down, ``1`` for up, and ``None`` for nearest, by default None

    Returns
    -------
    float
        rounded matrix values
    """        
    newval = np.zeros(np.shape(val))

    idx = (val != 0)
    power = np.floor(np.log10(abs(val[idx])))
    if toward == 0:
        newval[idx] = np.floor(val[idx] / 10**(power-order)) * 10**(power - order)
    elif toward == 1:
        newval[idx] = np.ceil(val[idx] / 10**(power-order)) * 10**(power - order)
    else:
        newval[idx] = np.round(val[idx] / 10**(power-order)) * 10**(power - order)

    return newval

def dynamic_format(value, threshold=1e4, order=4, toward=None):
    """Prepares number for display as *str*

    Formats a number of display as a *str*, typically in a *lineEdit* widget.

    Parameters
    ----------
    value : float
        number to round
    threshold : float, optional
        order of magnitude for determining display as floating point or expressing in engineering notation, by default 1e4
    order : int, optional
        number of orders to keep, by default 4
    toward : int, optional
        direction for rounding, ``0`` for down, ``1`` for up, and ``None`` for nearest, by default None

    Returns
    -------
    str
        number formatted as a string
    """        
    if toward is not None:
        value = oround(value, order=order, toward=toward)

    if abs(value) > threshold:
        return f'{{:.{order-1}e}}'.format(value)  # Scientific notation with order decimal places
    else:
        return f'{{:.{order}g}}'.format(value)


def symlog(x, linear_threshold=1.0):
    """Transforms data from linear to symlog space

    The advantage of a symlog transformation is that produces a log-like distribution for most
    values, but can be used on negative and zero values as well as positive values,
    preserves symmetry about 0.

    A symlog transformation is
     
    .. math::
        \operatorname{symlog} (x) = \operatorname{sign} (x) \log_{10} (1 + |x|/\lambda),

    where :math:`\lambda` is the linear threshold, region where the transformation is approximately linear.

    Parameters
    ----------
    x : numpy.ndarray
        the array of values to transform
    linear_threshold : float, optional
        linear threshold, often chosen as the mean, median, or standard deviation of the data
        or simply, by default 1.0

    Returns
    -------
    numpy.ndarray
        transformed values of x into symlog space
    """    
    x = np.asarray(x)
    return np.sign(x) * np.log10(1 + np.abs(x / linear_threshold))

def inv_symlog(x, linthresh=1.0):
    """Inverse symlog transform back to linear space.

    Converts a symlog array back to into linear space,

    .. math::
        \operatorname{symlog}^{-1}(x) = \operatorname{sign} (x) \lambda (10^{|x|} - 1)

    where :math:`\lambda` is the linear threshold used to control the region near zero where the
    data behave linearly.

    Parameters
    ----------
    x : numpy.ndarray
        the array of values to transform
    linear_threshold : float, optional
        linear threshold, often chosen as the mean, median, or standard deviation of the data
        or simply, by default 1.0

    Returns
    -------
    numpy.ndarray
        transformed values of x back to linear space
    :see also: symlog
    """    
    y = np.asarray(x)
    return np.sign(x) * linthresh * (10**np.abs(x) - 1)

def logit(x, eps=1e-8):
    """Applies the logit (log-odds) transformation to a numpy array.

    The logit function is defined as,

    .. math::
        \operatorname{logit}(x) = \log [x / (1 - x)]

    Values of `x` must be in the open interval (0, 1). To avoid division by zero or
    log of zero, values are clipped to [eps, 1 - eps].

    Parameters
    ----------
    x : numpy.ndarray
        Input array with values between 0 and 1.
    eps : float, optional
        Small value used to clip inputs away from 0 and 1 to prevent numerical issues.
        Default is 1e-8.

    Returns
    -------
    numpy.ndarray
        Transformed array with log-odds values.

    :see also: inv_logit
    """
    x = np.clip(x, eps, 1 - eps)
    return np.log(x / (1 - x))


def inv_logit(x):
    """Applies the inverse logit (sigmoid) transformation to a numpy array.

    The inverse logit transforms the data from normal space to sigmoid space and
    is defined as,

    .. math::
        \operatorname{logit}^{-1}(x) = [1 / (1 + e^{-x})]

    This maps real-valued inputs back to the (0, 1) interval.

    Parameters
    ----------
    x : numpy.ndarray
        Input array of real values.

    Returns
    -------
    numpy.ndarray
        Transformed array with values between 0 and 1.

    :see also: logit
    """
    return 1 / (1 + np.exp(-x))

def parse_isotope(field: str):
    """Converts an isotope into a symbol and mass.

    Separates an analyte field with a symbol-mass name to its separate parts.  For
    example, 27Al returns `Al` (symbol), `27` (mass).

    Parameters
    ----------
    field : str
        Field name to separate into symbol and mass.

    Returns
    -------
    symbol : str
        Returns the element symbol.
    mass : int
        Returns the isotope mass.
    """
    match = re.match(r"([A-Za-z]+)(\d*)", field)
    symbol = match.group(1) if match else field
    mass = int(match.group(2)) if match and match.group(2) else None

    return symbol, mass
