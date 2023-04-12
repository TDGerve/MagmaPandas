from .magmaFrame import MagmaFrame


class Clinopyroxene(MagmaFrame):
    @property
    def formula(self):
        """
        Docstrings
        """
        return self.mineral_formula(O=6)
