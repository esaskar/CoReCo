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

public class BondDistribution //implements Comparable<BondDistribution>
{
	public byte[] dist;
	
	public BondDistribution()
	{
	}
	
	public BondDistribution(Reaction r)
	{
		dist = new byte[Astar.bondstrsize];
		
		for (Bond b : r.getSubsBonds())
			dist[b.getBondStringIndex()]++;
		for (Bond b : r.getProdsBonds())
			dist[b.getBondStringIndex()]--;
	}
	
	public BondDistribution clone()
	{
		BondDistribution temp = new BondDistribution();
		temp.dist = new byte[Astar.bondstrsize];
		System.arraycopy(dist, 0, temp.dist, 0, Astar.bondstrsize);
		return temp;
	}
	
	public void update(AstarMapping m, Atom lhs, Atom rhs)
	{
		// update bond distribution when (lhs,rhs) pair is mapped
		// go through lhs neighbors and if any of the neighbors are included
		// in the map -> remove those from the distribution
		
		// update bond distribution
		for (Atom ne : lhs.getAtomNeighbors())
			if (ne.getOrder() < m.size())
				removeFromSub(ne.getBond(lhs));
		for (Atom ne : rhs.getAtomNeighbors())
			if (m.contains(ne))
				removeFromProd(ne.getBond(rhs));
	}
	
	public boolean equals(Object x)
	{
		return Arrays.equals(dist, ((BondDistribution)x).dist);
	}
	
	public int hashCode()
	{
		return Arrays.hashCode(dist);
	}
	
	public void removeFromSub(Bond b)
	{
		dist[b.getBondStringIndex()]--;
	}

	public void removeFromProd(Bond b)
	{
		dist[b.getBondStringIndex()]++;
	}
	
	public String toString()
	{
		return dist.toString();
	}

	public int bondCostEstimate()
	{
		int cost = 0;
		
		for (int i = 0; i < Astar.bondstrsize; i++)
			cost += Math.abs(dist[i]);
		
		return cost;
	}
}
