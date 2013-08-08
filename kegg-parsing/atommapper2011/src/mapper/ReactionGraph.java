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

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/*
 * ReactionGraph represents the reaction's atoms when mapped through a mapping
 *
 * The atoms and bonds of RG are distinct from 'Atom' and 'Bond' classes,
 *  but they contain references to them
 *  
 * RG is designed to be a separate entity from Reaction or Molecule objects
 * 
 * 
 */
enum Direction
{
	FORWARD,
	BACKWARD
}

public class ReactionGraph
{
	
	private Reaction reac;
	private Mapping m;
	
	private List<RGAtom> atoms = new ArrayList<RGAtom>();
	private List<RGBond> bonds = new ArrayList<RGBond>(); 
	
	
	// 'dir' is either +1 or -1 for forward and backward reactions, respectively
	public ReactionGraph(Reaction r, Mapping m)
	{
		this.reac = r;
		this.m = m;
		
		
		int[] index = new int[r.getSubsAtoms().size()];
		
		int id = 0;
		for (Atom a : r.getSubsAtoms())
		{
			index[a.getReactionalId()] = id;
			atoms.add(new RGAtom(id++, a.getSymbol(), this));
		}
		
		for (Atom a1 : r.getSubsAtoms())
		{
			for (Atom a2 : r.getSubsAtoms())
			{
				if (a1.getReactionalId() >= a2.getReactionalId())
					continue;
				
				Bond lhs = a1.getBond(a2);
				Bond rhs = m.getImage(a1).getBond( m.getImage(a2) );
				
				if (lhs != null || rhs != null)
				{
					RGBond b;
					int change = 0;
					
					if (lhs != null && rhs != null) // no-change
					{
						change = 0;
					}
					else if (lhs != null && rhs == null) // cleaved bond
					{
						change = -1;
					}
					else if (lhs == null && rhs != null) // new bond
					{
						change = +1;
					}
					
					RGAtom src = atoms.get(index[a1.getReactionalId()]);
					RGAtom tgt = atoms.get(index[a2.getReactionalId()]);
					int lhstype = lhs == null ? 0 : lhs.getType();
					int rhstype = rhs == null ? 0 : rhs.getType();
					
					b = new RGBond(src, tgt, change, lhstype, rhstype, this);
					bonds.add(b);
					
					src.addNeighbor(b, tgt);
					tgt.addNeighbor(b, src);
				}
			}
		}
	}
	
	
	
	// writes reaction graph as a mol file for both directions
	public void writeAsMol(String filename, Direction dir)
	{
		BufferedWriter out;

		String atomstr = rjust(Integer.toString(atoms.size()),3," ");
		String bondstr = rjust(Integer.toString(bonds.size()),3," ");
		
		try
		{
			out = new BufferedWriter(new FileWriter(filename));
			out.write(this.reac.getId() + "\n\n\n");
			out.write(atomstr + bondstr + "  0  0  0  0            999 V2000\n");

			// atoms
			for (RGAtom a : atoms)
				out.write("    0.0000    0.0000    0.0000 " + a.getSymbol() + "   0  0  0  0  0  0  0  0  0  0  0  0\n");
		
			// bonds
			for (RGBond b : bonds)
			{
				out.write(rjust(Integer.toString(b.getSource().getId()+1), 3, " "));
				out.write(rjust(Integer.toString(b.getTarget().getId()+1), 3, " "));
				if (dir == Direction.FORWARD)
				{
					out.write(rjust(Integer.toString(b.getChangetype()), 3, " "));
					out.write(rjust(Integer.toString(b.getOldtype()), 3, " "));
					out.write(rjust(Integer.toString(b.getNewtype()), 3, " "));
				}
				else if (dir == Direction.BACKWARD)
				{
					out.write(rjust(Integer.toString(-1*b.getChangetype()), 3, " "));
					out.write(rjust(Integer.toString(b.getNewtype()), 3, " "));
					out.write(rjust(Integer.toString(b.getOldtype()), 3, " "));
				}
				out.write( "  0  0\n");
			}
			
			out.write("M  END\n");
			out.close();
		}
		catch (IOException e)
		{
			System.out.println("X");
			System.out.println("ioerror: " + e.getMessage());
			return;
		}
	}
	
	public int size()
	{
		return atoms.size();
	}
	
	public List<RGAtom> getAtoms()
	{
		return atoms;
	}
	
	public List<RGBond> getBonds()
	{
		return bonds;
	}
	
	public String toString()
	{
		return "RG for " + reac;
	}
	
	public static String rjust(String s, int size, String pad)
	{
		while (s.length() < size)
			s = pad + s;
		return s;
	}
	
	public static String ljust(String s, int size, String pad)
	{
		while (s.length() < size)
			s = s + pad;
		return s;
	}
}



class RGAtom implements Comparable<RGAtom>
{
	private ReactionGraph parent;
	private int id;
	private String symbol;
	
	private Set<RGAtom> atomNeighbors = new HashSet<RGAtom>();
	private Map<RGAtom, RGBond> bondNeighbors = new HashMap<RGAtom, RGBond>();
	
	public RGAtom(int id, String symbol, ReactionGraph parent)
	{
		this.id = id;
		this.symbol = symbol;
		this.parent = parent;
	}
	
	public String getSymbol()
	{
		return symbol;
	}
	
	public int getId()
	{
		return id;
	}
	
	public ReactionGraph getParent()
	{
		return parent;
	}
	
	public Collection<RGAtom> getAtomNeighbors()
	{
		return atomNeighbors;
	}

	public boolean isNeighbor(RGAtom other)
	{
		return getBond(other) != null;
	}
	
	public RGBond getBond(RGAtom other)
	{
		return bondNeighbors.get(other);
	}
	
	public void addNeighbor(RGBond b, RGAtom a)
	{
		atomNeighbors.add(a);
		bondNeighbors.put(a, b);
	}
		
	public String toString()
	{
		return id + "/" + symbol;
	}

	public boolean equals(Object o)
	{
		RGAtom other = (RGAtom)o;
		return id == other.getId() && parent == other.getParent();
	}
	
	public int hashCode()
	{
		return id;
	}

	public int compareTo(RGAtom o)
	{
		if (id < o.getId())
			return -1;
		else if (id > o.getId())
			return 1;
		return 0;
	}
}



class RGBond
{
	private RGAtom source;
	private RGAtom target;
	private int changetype; // the change pattern, 0 (no change), +1 (new bond), -1 (cleaved bond)
	private int oldtype; // the bond's old type, e.g. 2 for double bond
	private int newtype; // the bond's new type, e.g. 0 for no bond (when cleavage happens)
	private ReactionGraph parent;
	
	public RGBond(RGAtom source, RGAtom target, int change, int oldtype, int newtype, ReactionGraph parent)
	{
		this.source = source;
		this.target = target;
		this.changetype = change;
		this.oldtype = oldtype;
		this.newtype = newtype;
		this.parent = parent;
	}
	
	public RGAtom getSource()
	{
		return source;
	}

	public RGAtom getTarget()
	{
		return target;
	}

	public int getChangetype()
	{
		return changetype;
	}

	public int getOldtype()
	{
		return oldtype;
	}

	public int getNewtype()
	{
		return newtype;
	}
	
	public ReactionGraph getParent()
	{
		return parent;
	}
	
	public String toString()
	{
		return source.getId() + "<->" + target.getId() + "(" + changetype + ":" + oldtype + "->" + newtype + ")";
	}
}


