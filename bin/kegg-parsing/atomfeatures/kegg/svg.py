#!/usr/bin/python

# some functions to support molecular drawings using svg

class Svg:
	def __init__(self):
		self.lines = []
		self.s = ""
		self.width = 0
		self.height = 0
		self.points = []
	
	def header(self):
		lines.insert(0, "<?xml version=\"1.0\"?>")
		lines.insert(1, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">")
		lines.insert(2, "")
		lines.insert(3, "<svg fill-opacity=\"1\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" color-rendering=\"auto\" color-interpolation=\"auto\" stroke=\"black\" text-rendering=\"auto\" stroke-linecap=\"square\" width=\"%d\" stroke-miterlimit=\"10\" stroke-opacity=\"1\" shape-rendering=\"auto\" fill=\"black\" stroke-dasharray=\"none\" font-weight=\"normal\" stroke-width=\"1\" height=\"%d\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"&apos;Dialog&apos;\" font-style=\"normal\" stroke-linejoin=\"miter\" font-size=\"12\" stroke-dashoffset=\"0\" image-rendering=\"auto\">" % (self.width, self.height))
		lines.insert(4, "<g>")
	
	def footer(self):
		lines.append("</g>")
		lines.append("</svg>")
	
	# public method
	def add_symbol(self, text, x, y, color="black"):
		lines.append("<rect x=\"%d\" y=\"%d\" fill=\"white\" width=\"15\" height=\"16\" stroke=\"none\" />" % (x-7,y-8))
		lines.append("<text fill=\"%s\" x=\"%d\" y=\"%d\" stroke=\"none\">%s</text>" % (color,x-5,y+5,text))
		self.points.append((x,y))
	
	# public method
	def add_line(self, x1, y1, x2, y2):
		lines.append("<line x1=\"%f\" fill=\"none\" y1=\"%f\" x2=\"%f\" y2=\"%f\" />" % (x1,y1,x2,y2))
	
	def count_size(self):
		maxx = max([p[0] for p in self.points])
		maxy = max([p[1] for p in self.points])
		
		return maxx,maxy
	
	def finalize(self):
		self.width, self.height = self.count_size()
		self.header()
		self.footer()
		
		self.s = ""
		prefix = 0
		for line in lines:
			self.s += (" "*prefix) + line + "\n"
			
			if "<g>" in line:
				prefix += 1
			elif "</g>" in line:
				prefix -= 1
	
	# public method
	def write(self, fname):
		if not self.s:
			self.finalize()
		
		f = open(fname, "w")
		f.write(self.s)
		f.close()
	
	def __str__(self):
		if not self.s:
			self.finalize()
		return self.s






