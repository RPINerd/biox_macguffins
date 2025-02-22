"""All the classes that are used across the biox_macguffins codebase"""


class Primer:

    """A class to represent a primer"""

    def __init__(self, chr: str | int, num: int, start: int, end: int) -> None:
        """Initializes a primer object"""
        self.chr = chr
        self.number = num
        self.start = start
        self.end = end

    def __str__(self) -> str:
        """Returns a string representation of the primer"""
        return f"Primer {self.number} on chromosome {self.chr} from {self.start} to {self.end}"

    def __repr__(self) -> str:
        """Returns a string representation of the primer"""
        return f"Primer{self.number}, {self.chr}:{self.start}-{self.end}"
