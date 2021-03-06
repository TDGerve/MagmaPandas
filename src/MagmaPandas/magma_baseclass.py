from typing import List
import pandas as pd
import elements as e


def read_file(file: str, *args, index_col: List[str], keep_columns: List[str] = [], units=None, type=None, **kwargs):
    """
    Docstring
    """
    allowed_units = ['wt. %', 'mol fraction']
    allowed_types = ['oxide', 'cation']
    if not units in allowed_units:
        raise ValueError(f'units: "{units}" not recognised, please choose from: {allowed_units}')
    if type not in allowed_types:
        raise ValueError(f'type: "{type}" not recognised, please choose from: {allowed_types}')

    df = pd.read_csv(file, index_col=index_col, **kwargs)

    return MagmaBase(df, *args, keep_columns=keep_columns, calculate_total=True, units=units, datatype=type, **kwargs)


class MagmaBase(pd.DataFrame):

    _metadata = ["_weights", "_no_data", "_units", "_datatype"]

    def __init__(
        self,
        df: pd.DataFrame,
        *args,
        keep_columns: List[str]=[],
        calculate_total=False,
        units=None,
        datatype=None,
        **kwargs,
    ):
        # A pandas series with the masses of all oxides and elements in the dataframe
        self._weights = pd.Series(name="weight", dtype="float32")
        # A list with the names of all columns that do not contain chemical data
        self._no_data = keep_columns
        self._units = units
        self._datatype = datatype

        if isinstance(df, pd.DataFrame):
            
            delete_columns = set()
            keep_columns = set(keep_columns)

            for col in df.columns:
                try:
                    self._weights[col] = e.calculate_weight(col)
                except:
                    delete_columns.add(col)

            # Drop all columns without chemical data, unless explicitely specified otherwise in 'keep_columns'
            df = df.drop(delete_columns.difference(keep_columns), axis=1)
            df[self.elements] = df[self.elements].astype('float32')
            # Calculate total concentrattions
            if calculate_total:
                df["total"] = df[self.elements].sum(axis=1)

        super().__init__(df, *args, **kwargs)


    @property
    def _constructor(self):
        """This is the key to letting Pandas know how to keep
        derivatives of `MagmaBase` the same type as yours.  It should
        be enough to return the name of the Class.  However, in
        some cases, `__finalize__` is not called and `new attributes` are
        not carried over.  We can fix that by constructing a callable
        that makes sure to call `__finalize__` every time."""

        def _c(*args, **kwargs):
            return MagmaBase(*args, **kwargs).__finalize__(self)
        return _c


    @property
    def units(self):
        return f'{self._datatype} {self._units}'

    @units.setter
    def units(self, value):
        print('units are read only')


    @property
    def weights(self):
        """
        Docstrings
        """        
        return self._weights


    @property
    def elements(self):
        """
        Docstrings
        """        
        return self._weights.index.values


    @property
    def moles(self):
        """
        Docstrings
        """
        if self._datatype != 'oxide':
            raise TypeError(f'{self} is not in oxides')
        if self._units != 'wt. %':
            raise TypeError(f'{self} is not in wt. %')
        moles = self.copy()[self.elements]
        # Calculate moles
        moles = moles.div(moles._weights)
        # Normalise to 1
        total = moles.sum(axis=1)
        moles = moles.div(total, axis=0)
        moles["total"] = moles.sum(axis=1)
        # Add back columns without chemical data
        moles[self._no_data] = self[self._no_data]
        moles._units = 'mol fraction'

        return moles


    @property
    def cations(self):
        """
        Docstrings
        """
        if self._datatype != 'oxide':
            raise TypeError(f'{self} is not in oxides')
        # Calculate oxide moles
        moles = self.moles[self.elements]
        # Calculate cation moles
        cation_numbers = e.cation_numbers(moles.elements)
        cations = moles.mul(cation_numbers)
        # Rename columns to cations
        cations.columns = e.cation_names(cations.elements)
        # Normalise to 1
        total = cations.sum(axis=1)
        cations = cations.div(total, axis=0)
        cations['total'] = cations.sum(axis=1)
        # Add back columns without chemical data
        cations[self._no_data] = self[self._no_data]
        # Set datatype and elements
        cations._datatype = 'cation'
        cations.recalculate()

        return cations


    def mineral_formula(self, O: int=None):
        """
        Docstrings
        """
        # Calculate cation fractions
        cations = self.cations
        cations = cations[cations.elements]
        # Calculate oxygens per cation
        oxygen_numbers = e.oxygen_numbers(self.elements) / e.cation_numbers(self.elements)
        oxygen_numbers.index = cations.elelements
        # Normalise to oxygen
        oxygen_total = cations.mul(oxygen_numbers).sum(axis=1)
        oxygen_factor = O / oxygen_total
        cations = cations.mul(oxygen_factor, axis=0)
        cations['O'] = O

        return cations

    
    def recalculate(self):
        """
        Docstrings
        """
        weights = pd.Series(name="weight", dtype="float32")

        for col in self.columns:
                try:
                    weights[col] = e.calculate_weight(col)
                except:
                    pass
                 
        self._weights = weights

        self['total'] = self[self.elements].sum(axis=1)