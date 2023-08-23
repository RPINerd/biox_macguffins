"""
    All the classes that are used across the biox_macguffins codebase
"""


class Primer:
    def __init__(self, chr, num, start, end):
        self.chr = chr
        self.number = num
        self.start = start
        self.end = end

    def chr(self):
        return self.chr

    def start(self):
        return self.start

    def end(self):
        return self.end

    def number(self):
        return self.number
