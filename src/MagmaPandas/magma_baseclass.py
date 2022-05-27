from typing import List
import pandas as pd
import elements as e


def read_file(file: str, index_col: List[str], keep_columns: List[str] = [], **kwargs):
    """
    Docstring
    """

    df = pd.read_csv(file, index_col=index_col, **kwargs)

    return MagmaBase(df, keep_columns=keep_columns, calculate_total=True, **kwargs)


class MagmaBase(pd.DataFrame):

    _metadata = ["_weights", "_no_data", "_units"]

    @property
    def _constructor(self):
        """This is the key to letting Pandas know how to keep
        derivative `MagmaBase` the same type as yours.  It should
        be enough to return the name of the Class.  However, in
        some cases, `__finalize__` is not called and `new attributes` are
        not carried over.  We can fix that by constructing a callable
        that makes sure to call `__finalize__` every time."""

        def _c(*args, **kwargs):
            return MagmaBase(*args, **kwargs).__finalize__(self)
        return _c

    def __init__(
        self,
        df: pd.DataFrame,
        *args,
        keep_columns: List[str]=[],
        calculate_total=False,
        units='wt. %',
        **kwargs,
    ):
        # A pandas series with the masses of all oxides and elements in the dataframe
        self._weights = pd.Series(name="weight", dtype="float32")
        # A list with the names of all columns that do not contain chemical data
        self._no_data = keep_columns
        self._units = units


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
    def units(self):
        return self._units

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
        moles = self.copy()[self.elements]
        # Calculate moles
        moles = moles.div(moles._weights)
        # Normalise to 1
        total = moles.sum(axis=1)
        moles = moles.div(total, axis=0)
        moles["total"] = moles.sum(axis=1)
        # Add back columns without chemical data
        moles[self._no_data] = self[self._no_data]
        moles._units = 'moles'
        return moles


    @property
    def cations(self):
        """
        Docstrings
        """
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
        cations._units = 'cations'
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