class Linebuffer:
    def __init__(self, f):
        self.f = f
        self.buf = None
    def getLine(self):
        if self.buf == None:
            self.advance()
        return self.buf
    def advance(self):
        self.buf = self.f.readline()
