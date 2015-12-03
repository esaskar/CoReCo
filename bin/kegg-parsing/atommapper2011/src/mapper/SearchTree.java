//
//    Copyright 2011 Markus Heinonen 
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
package mapper;

import java.io.*;

public class SearchTree
{
	private StringBuffer treebuf = new StringBuffer();
	
	public SearchTree(Atom[] subatoms)
	{
		insertHeader();
		
		treebuf.append("\"left:root\" -> " + "\"left:" + Integer.toString(subatoms[0].getOrder()) + ":" + subatoms[0].getSymbol() + "\";\n");
		
		for (int i = 0; i < subatoms.length - 1; i++)
			treebuf.append("\"left:" + Integer.toString(subatoms[i].getOrder()) + ":" + subatoms[i].getSymbol() + "\" -> \"left:" + Integer.toString(subatoms[i+1].getOrder()) + ":" + subatoms[i+1].getSymbol() + "\";\n");
	}
	
	private void insertHeader()
	{
		treebuf.append("digraph G {\n");
		treebuf.append("node [style=dotted];\n");
	}
	
	private void insertFooter()
	{
		treebuf.append("\n}");
	}
	private void removeFooter()
	{
		treebuf = treebuf.delete(treebuf.length()-1, treebuf.length());
	}
	
	public void addNode(AstarMapping prev, AstarMapping m)
	{
		String prevrow;
		String prevsymbol;
		String prevcost;
		
		if (prev.size() == 0)
		{
			prevrow = "root";
			prevsymbol = "";
			prevcost = "";
		}
		else
		{
			prevrow = Integer.toString( prev.getMap()[ prev.size()-1 ] );
			prevsymbol = Astar.subatoms[ prev.getMap()[ prev.size()-1 ] ].getSymbol();
			prevcost = "(" + prev.getGCost() + "+" + prev.getHCost() + ")";
		}

		String rownumber = Integer.toString( m.getMap()[ m.size()-1 ] ); 
		String symbol = Astar.subatoms[ m.getMap()[ m.size()-1 ] ].getSymbol();
		String cost = "(" + m.getGCost() + "+" + m.getHCost() + ")";
//		String featcost = "(" + (m.getFeatDiff()*100) + ")";
		
		if (prevrow.equals("root"))
			treebuf.append("\"root\" -> \"" + rownumber + ":" + symbol + ":" + m.size() + cost + "\";\n");
		else
			treebuf.append("\"" + prevrow + ":" + prevsymbol + ":" + prev.size() + prevcost + "\" -> \"" + rownumber + ":" + symbol + ":" + m.size() + cost + "\";\n");
	}
	
	public void chooseNode(AstarMapping m)
	{
		String rownumber = Integer.toString( m.getMap()[ m.size()-1 ] ); 
		String symbol = Astar.subatoms[ m.getMap()[ m.size()-1 ] ].getSymbol();
		String cost = "(" + m.getGCost() + "+" + m.getHCost() + ")";
		
		treebuf.append("\"" + rownumber + ":" + symbol + ":" + m.size() + cost + "\" [style=solid]; \n");
	}
	
	private void writeTree()
	{
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter("searchtree.dot"));
			out.write(treebuf.toString());
			out.close();
		} catch (IOException e)
		{
			System.out.println(e.getMessage());
		}
	}
	
	public void drawTree()
	{
		insertFooter();
		writeTree();
		
		try {
			Runtime.getRuntime().exec("dot -Tsvg -o searchtree.svg searchtree.dot");
		} catch (IOException e)
		{
			System.out.println(e.getMessage());
		}
		
		removeFooter();
	}
}

