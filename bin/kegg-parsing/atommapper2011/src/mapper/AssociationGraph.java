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


public class AssociationGraph
{
	@SuppressWarnings("unused")
	private Reaction reac;
	
	private Collection<AGNode> nodes;
	private Collection<AGEdge> edges;
	
	public AssociationGraph(Reaction r)
	{
		reac = r;
		
		r.computeDistances();
		
		this.reac = r;
		this.nodes = new HashSet<AGNode>();
		this.edges = new HashSet<AGEdge>();
		
		Collection<Atom> sub = r.getSubsAtoms();
		Collection<Atom> prod = r.getProdsAtoms();
		
		for (Atom lhs : sub)
			for (Atom rhs : prod)
				if (lhs.getSymbol().equals(rhs.getSymbol()))
					nodes.add(new AGNode(lhs, rhs));

		for (AGNode agn1 : nodes)
		{
			for (AGNode agn2 : nodes)
			{
				if (agn1 == agn2 || agn1.getLHS() == agn2.getLHS() || agn1.getRHS() == agn2.getRHS())
					continue;
				
				if (r.distance(agn1.getLHS(), agn2.getLHS()) < 0 || r.distance(agn1.getRHS(), agn2.getRHS()) < 0)
					continue;
				
				boolean agn1neigh = agn1.getLHS().isNeighbor(agn2.getLHS());
				boolean agn2neigh = agn1.getRHS().isNeighbor(agn2.getRHS());
				
				if (agn1neigh == agn2neigh)
				{
					edges.add(new AGEdge(agn1, agn2));
					
					agn1.addNeighbor(agn2);
					agn2.addNeighbor(agn1);
				}
			}
		}
	}
	
	public Collection<AGNode> getNodes()
	{
		return nodes;
	}
}


