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

public class Bond
{
	// private int id;
	private Atom a1;
	private Atom a2;
	private int type;
	private Molecule parent;
	private int bondStringIndex;
	private String bondString;
	
	public static final int SINGLE_BOND = 1;
	public static final int DOUBLE_BOND = 2;
	public static final int TRIPLE_BOND = 3;

	public Bond(Atom a1, Atom a2, int type, boolean update_nbs)
	{
		// a1 << a2, based on the rownumber
		if (a1.getMolecularId() > a2.getMolecularId())
		{
			Atom apu = a2;
			a2 = a1;
			a1 = apu;
		}

		this.a1 = a1;
		this.a2 = a2;
		this.type = type;
		
		if (a1.getSymbol().compareTo(a2.getSymbol()) <= 0)
		{
			if (GlobalOptions.bc)
				bondString = a1.getSymbol() + type + a2.getSymbol();
			else
				bondString = a1.getSymbol() + a2.getSymbol();
		}
		else
		{
			if (GlobalOptions.bc)
				bondString = a2.getSymbol() + type + a1.getSymbol();
			else
				bondString = a2.getSymbol() + a1.getSymbol();
		}
		
		if (update_nbs)
		{
			this.a1.addNeighbor(this, a2);
			this.a2.addNeighbor(this, a1);
		}
	}

	public String getBondString()
	{
		return bondString;
	}
	
	public int getBondStringIndex()
	{
		return bondStringIndex;
	}
	
	public void setBondStringIndex(int x)
	{
		bondStringIndex = x;
	}
	
	public int getType()
	{
		return type;
	}

	public Atom getA1()
	{
		return a1;
	}

	public Atom getA2()
	{
		return a2;
	}

	public void setParent(Molecule m)
	{
		this.parent = m;
	}

	public Molecule getParent()
	{
		return this.parent;
	}

	public boolean equals(Object bo)
	{
		Bond b = (Bond)bo;
		Atom at1 = b.getA1();
		Atom at2 = b.getA2();

		if (this.a1.equals(at1) && this.a2.equals(at2)
				&& this.type == b.getType())
			return true;
		if (this.a1.equals(at2) && this.a2.equals(at1)
				&& this.type == b.getType())
			return true;

		return false;
	}
	
	public String toString()
	{
		String ret = "" + this.a1.getSymbol() + "(" + this.a1.getMolecularId()
				+ "/" + this.a1.getParent().getMolNum() + ")";
		ret += this.type + this.a2.getSymbol() + "(" + this.a2.getMolecularId()
				+ "/" + this.a2.getParent().getMolNum() + ")";
		return ret;
	}
	
	public String stringRepr()
	{
		if (a1.getSymbol().compareTo(a2.getSymbol()) <= 0)
			return new String(a1.getSymbol() + this.type + a2.getSymbol());
		else
			return new String(a2.getSymbol() + this.type + a1.getSymbol());		
	}

//	public static Collection<Bond> setMinus(Collection<Bond> a, Collection<Bond> b)
//	{
//		Collection<Bond> x = new HashSet<Bond>();
//		x.addAll(a);
//		x.removeAll(b);
//		return x;
//	}
}
