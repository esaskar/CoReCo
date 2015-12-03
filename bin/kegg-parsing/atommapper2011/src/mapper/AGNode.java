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

public class AGNode implements Comparable<AGNode>
{
	private Atom lhs;
	private Atom rhs;
	private Collection<AGNode> neighs;
	
	public AGNode(Atom a, Atom b)
	{
		neighs = new HashSet<AGNode>();
		lhs = a;
		rhs = b;
	}
	
	public Atom getLHS()
	{
		return lhs;
	}
	public Atom getRHS()
	{
		return rhs;
	}
	
	public boolean isNeighbor(AGNode other)
	{
		return neighs.contains(other);
	}
	
	public void addNeighbor(AGNode other)
	{
		neighs.add(other);
	}
	
	public boolean neighborPair()
	{
		return lhs.isNeighbor(rhs);
	}
	
	public String toString()
	{
		return lhs.toString() + " <=> " + rhs.toString();
	}

	public int compareTo(AGNode o)
	{
		if (lhs.getReactionalId() < o.getLHS().getReactionalId())
			return -1;
		else if (lhs.getReactionalId() == o.getLHS().getReactionalId() && rhs.getReactionalId() < o.getRHS().getReactionalId())
			return -1;
		else if (lhs.getReactionalId() == o.getLHS().getReactionalId() && rhs.getReactionalId() == o.getRHS().getReactionalId())
			return 0;
		
		return 1;
	}
	
	public int hashCode()
	{
		// returns something 
		return lhs.hashCode() * rhs.hashCode();
	}
	
	public boolean equals(Object other)
	{
		AGNode ot = (AGNode)other;
		
		return lhs == ot.getLHS() && rhs == ot.getRHS();
	}
	
}
