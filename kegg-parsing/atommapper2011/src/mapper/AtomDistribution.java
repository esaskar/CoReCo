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

import java.util.Arrays;


// Represents Distribution of atomNeighStrings on both sides of the reaction mapping
// e.g. "C|CCO:+2", "O|CC:-1", ...
public class AtomDistribution
{
	private byte[] dist2;
	
	public AtomDistribution()
	{
		dist2 = new byte[Astar.atomstrsize];
	}

	public AtomDistribution(Reaction r)
	{
		dist2 = new byte[Astar.atomstrsize];
		
		for (Atom a : r.getSubsAtoms())
			dist2[a.getNeighStringInt()]++;
		for (Atom a : r.getProdsAtoms())
			dist2[a.getNeighStringInt()]--;
	}
	
	public int atomCostEstimate()
	{
		int cost = 0;
		
		for (int i = 0; i < Astar.atomstrsize; i++)
			if (dist2[i] > 0)
				cost += dist2[i];
				
		// should work!
		return (cost + 1) / 2;
	}

	public AtomDistribution clone()
	{
		AtomDistribution temp = new AtomDistribution();
		temp.dist2 = new byte[Astar.atomstrsize];
		System.arraycopy(dist2, 0, temp.dist2, 0, Astar.atomstrsize);
		return temp;
	}
	
	public void update(Atom lhs, Atom rhs)
	{
		int lhss = lhs.getNeighStringInt();
		int rhss = rhs.getNeighStringInt();
		
		dist2[lhss]--; // we are removing atoms from unmapped region -> reverse direction
		dist2[rhss]++;
	}
	
	public String toString()
	{
		return dist2.toString();
	}
	
	public boolean equals(Object x)
	{
		return Arrays.equals(dist2, ((AtomDistribution)x).dist2);
	}
	
	public int hashCode()
	{
		return Arrays.hashCode(dist2);
	}
}



