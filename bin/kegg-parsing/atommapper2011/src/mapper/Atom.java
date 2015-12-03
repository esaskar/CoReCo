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

import java.util.*;

public class Atom implements Comparable<Atom>
{
	private Molecule parent;
	private String symbol = "";
	// unique id of single molecule atoms
	private int id_of_molecule; // 0...n, where n=atomcount of molecule
	// Unique id over reactions
	private int id_of_reaction; // 0...m, where m=atomcount of reaction
	// unique id over rhs or lhs
	private int order = -1;
	// atom neighborhood string, ex C|COO for carbon with two O's and one C neighbor
	private String atomNeighString = null;
	private int atomNeighStringIndex = -1;
	
	private Set<Atom> atomNeighbors = new HashSet<Atom>();
	private Map<Atom, Bond> bondNeighbors = new HashMap<Atom, Bond>();

	// product side atoms are organized into bins
	// each atom's location is prodatoms[bin_id][bin_offset]
	public  int bin_id;
	public  int bin_offset;
	
	
	// create a node without any information,
	// don't increment ids
	public Atom()
	{
		this.id_of_reaction = -1;
		this.id_of_molecule = -1;
	}

	public Atom(int id, String c)
	{
		this.symbol = c;
		this.id_of_molecule = id;
		this.id_of_reaction = -1;
	}

	public Collection<Atom> getAtomNeighbors()
	{
		return this.atomNeighbors;
	}
	
	public Collection<Bond> getBondNeighbors()
	{
		return this.bondNeighbors.values();
	}
	
	public String getSymbol()
	{
		return this.symbol;
	}

	private void setSymbol(String c)
	{
		this.symbol = c;
	}

	public int getMolecularId()
	{
		return id_of_molecule;
	}
	
	public void setMolecularId(int i)
	{
		id_of_molecule = i;
	}
	
	public Collection<Atom> getNeighbours()
	{
		return this.atomNeighbors;
	}

	public boolean isNeighbor(Atom other)
	{
		return getBond(other) != null;
	}
	
	public Bond getBond(Atom other)
	{
		return bondNeighbors.get(other);
	}

	private void setNeighbors(Map<Atom, Bond> c)
	{
		for (Atom a : c.keySet())
			addNeighbor(c.get(a), a);
	}

	public void addNeighbor(Bond b, Atom a)
	{
//		if (!atomNeighbors.contains(a))
			atomNeighbors.add(a);

//		if (!bondNeighbors.contains(b))
			bondNeighbors.put(a, b);
	}
	
	public void removeNeighbor(Atom a)
	{
		this.atomNeighbors.remove(a);
		this.bondNeighbors.remove(a);
	}

	public int getDegree()
	{
		return atomNeighbors.size();
	}
	
	public void setParent(Molecule m)
	{
		this.parent = m;
	}

	public Molecule getParent()
	{
		return this.parent;
	}

	public int getReactionalId()
	{
		return this.id_of_reaction;
	}
	
	public void setReactionalId(int value)
	{
		id_of_reaction = value;
	}

	public String toString()
	{
		return this.symbol + "(" + this.id_of_molecule + "/" + this.parent.getMolNum() + ")";
	}

	// requires id's to be unique!
	public boolean equals(Object at)
	{
		return this.id_of_reaction == ((Atom)at).getReactionalId() && this.parent == ((Atom)at).getParent();
	}

	public int hashCode()
	{
		return this.id_of_reaction;
	}
	
	public Atom clone()
	{
		Atom a = new Atom();
		a.setParent(this.parent);
		a.setMolecularId(this.id_of_molecule);
		a.setReactionalId(this.id_of_reaction);
		a.setOrder(this.order);
		a.setSymbol(this.symbol);
		a.setNeighbors(this.bondNeighbors);
		return a;
	}

	public int getOrder()
	{
		return order;
	}
	
	public void setOrder(int value)
	{
		this.order = value;
	}
	
	public String getNeighString()
	{
		if (atomNeighString != null)
			return atomNeighString;
		
		String astr = symbol + "|";
		String nstr = "";
		Collection<Atom> neighs = getNeighbours();
		for (Atom ne : neighs)
			nstr += ne.getSymbol();
		
		char[] chars = nstr.toCharArray();
		Arrays.sort(chars);
		astr += new String(chars);
		
		atomNeighString = astr;

		return atomNeighString;
	}
	
	// each neighstring (e.g. "C|COO") is indexed so that the set of neighstring
	// can be put into a regular array.
	// returns this "hash number"
	public int getNeighStringInt()
	{
		return atomNeighStringIndex;
	}
	
	public void setNeighStringInt(int x)
	{
		atomNeighStringIndex = x;
	}
	
	public int compareTo(Atom a)
	{
		if (this.order < a.getOrder())
			return -1;
		else if (this.order > a.getOrder())
			return 1;
		return 0;
	}

//	// assume a >= b
//	public static Collection<Atom> setMinus(Collection<Atom> a, Collection<Atom> b)
//	{
//		Collection<Atom> ret = new HashSet<Atom>();
//		for (Atom at : a)
//		{
//			if (!b.contains(at))
//				ret.add(at);
//		}
//		return ret;
//	}
//
//	public static boolean setEquals(Collection<Atom> A, Collection<Atom> B)
//	{
//		for (Atom a : A)
//			if (!B.contains(a))
//				return false;
//		for (Atom b : B)
//			if (!A.contains(b))
//				return false;
//		return true;
//	}
}
