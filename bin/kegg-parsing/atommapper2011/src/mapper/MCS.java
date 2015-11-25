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

public class MCS
{
	private static int round = 0;
	private static LinkedList<AGNode> current = new LinkedList<AGNode>();
	private static Collection<AGNode[]> cliques = new HashSet<AGNode[]>();
	
	public static List<Mapping> mcs(Reaction r)
	{
		ArrayList<Mapping> maps = new ArrayList<Mapping>();
		Mapping map = new HashMapping();
				
		// iteratively chop MCS-block away from reaction
		while (!r.getSubsAtoms().isEmpty())
		{
			// build the association graph of current reaction graph,
			//  find maximum clique, turn it into a maximum common subgraph (ie. mapping),
			//  and remove those parts from reaction graph
			AssociationGraph ag = new AssociationGraph(r);
			Collection<AGNode[]> maxcliqs = findMCL(ag);
			HashMapping cliqmap = MCLtoMCS(maxcliqs);
			r.remove(cliqmap);
			
			// add cliqmap to mapping
			for (Atom a : cliqmap.getDomain())
				map.extend(a, cliqmap.getImage(a));
		}
		
		r.restore();
		
		map.setGCost(map.costFunction());
//		map.setFCost();
		
		maps.add(map);
		
		return maps;
	}
	
	private static Collection<AGNode[]> findMCL(AssociationGraph ag)
	{
		cliques.clear();
//		current.clear();
		
		findGreedyCliques(ag);
		findCliques(ag);
		
		return cliques;
	}

	// lots of same-size mcs's, all have zero-cost
	// pick first one which is connected (due to some bug sometimes they are disjoint)
	private static HashMapping MCLtoMCS(Collection<AGNode[]> cliques)
	{
		HashMapping best = null;
		
		for (AGNode[] cliq : cliques)
		{
			HashMapping m = new HashMapping();
			
			for (AGNode a : cliq)
				m.extend(a.getLHS(), a.getRHS());
			
			if (m.connected())
				best = m;
		}
		
		if (best == null)
			return getMaxConnectedMapping(cliques);
		else
			return best;
	}
	
	// all cliques are disjoint
	private static HashMapping getMaxConnectedMapping(Collection<AGNode[]> cliques)
	{
		HashMapping best = null;
		System.out.println("Only disjoint mcs's found, extracting largest piece from them");
		
		for (AGNode[] cliq : cliques)
		{
			// create all connected mappings of a clique, try to place largest to 'best'
			HashMapping m = new HashMapping();
			
			for (AGNode a : cliq)
				m.extend(a.getLHS(), a.getRHS());
			
			Collection<HashMapping> pieces = m.pieces();
			
			for (HashMapping p : pieces)
				if (best == null || p.size() > best.size())
					best = p;
		}
		
		return best;
	}
	
	
	// finds fast a set of candidate cliques by starting the search from each node
	private static void findGreedyCliques(AssociationGraph ag)
	{
//		LinkedList<AGNode> current = new LinkedList<AGNode>();
		current.clear();
		HashSet<AGNode> border = new HashSet<AGNode>(ag.getNodes()); // full
		HashSet<AGNode> closed = new HashSet<AGNode>(); // empty
		
		for (AGNode first : border)
		{
			closed.clear();
			current.clear();
			current.add(first);
			
			HashSet<AGNode> new_border = new HashSet<AGNode>();
			for (AGNode b : ag.getNodes())
				if (first.isNeighbor(b))// && first.getLHS().isNeighbor(b.getLHS()) && first.getRHS().isNeighbor(b.getRHS()))// && first.) // connected in original graphs
					new_border.add(b);
			
			round = 0;
			// only run until border exhausted
			extend(new_border, closed, new_border.size()+2);
		}
	}	
	
	private static void findCliques(AssociationGraph ag)
	{
		// uses three lists
		// - 'current' is the current clique-in-construction
		// - 'border' are nodes completely connected to the 'compsub'
		// - 'closed' are nodes already processed which lead to valid extensions

		
//		LinkedList<AGNode> current = new LinkedList<AGNode>();
		current.clear(); 
		HashSet<AGNode> border = new HashSet<AGNode>(ag.getNodes()); // full
		HashSet<AGNode> closed = new HashSet<AGNode>(); // empty
		
		round = 0;
		extend(border, closed, 10000);
	}
	
	private static void extend(HashSet<AGNode> border, HashSet<AGNode> closed, int limit)
	{
		round++;
		
		if (round > limit)
			return;
			
		// (1) select a candidate
		// (2) add it to compsub
		// (3) create new candidates and noway
		// (4) call extend on those
		// (5) move selected candidate from compsub to noway

		Iterator<AGNode> it = border.iterator();
		boolean noneConnected = true;
		
		while (!border.isEmpty() && it.hasNext() && round <= limit)
		{
			if (!it.hasNext())
				System.out.println("something fishy is going on");
				
			AGNode cand = it.next();
//			AGNode cand = border.iterator().next();
			
			// check whether current is connected, if not, restore previous current and continue
			if (connected(current, cand) == false)
				continue;
			
			noneConnected = false;
			it.remove(); // remove cand from border cleanly
//			border.remove(cand);
			current.addLast(cand);
			
			HashSet<AGNode> new_closed = new HashSet<AGNode>();
			for (AGNode a : closed)
				if (cand.isNeighbor(a))
					new_closed.add(a);
			
			HashSet<AGNode> new_border = new HashSet<AGNode>();
			for (AGNode a : border)
				if (cand.isNeighbor(a))
					new_border.add(a);
			
			extend(new_border, new_closed, limit);
			
			current.removeLast(); // last
			closed.add(cand);
		}
	
		// found clique if nothing to add
		if ( (border.isEmpty() || noneConnected) && closed.isEmpty())
		{
			// halutaan maximaalinen klikki, sen jälkeen se joka minimoi costfunktion...
			
			AGNode[] clique = new AGNode[current.size()];
			current.toArray(clique);
			Arrays.sort(clique);
			
			if (cliques.isEmpty() || clique.length > cliques.iterator().next().length)
			{
				cliques.clear();
				cliques.add(clique);
			}
			else if (clique.length == cliques.iterator().next().length)
			{
				boolean contains = false;
				// test if result already there...
				for (AGNode[] res : cliques)
				{
					contains = true;
					
					if (clique.length != res.length)
					{
						contains = false;
						continue;
					}
					
					for (int i = 0; i < res.length; i++)
					{
						if (res[i] != clique[i])
						{
							contains = false;
							break;
						}
					}
					
					if (contains)
						break;
				}
	
				if (!contains)
					cliques.add(clique);
			}
		}				
		
	}
	
	private static boolean connected(Collection<AGNode> agnodes, AGNode candidate)
	{
		// checks whether 'candidate' AGNode is connected from left and right to corresponding graphs...
		// simple
		
		if (agnodes.size() == 0)
			return true;		
		
		boolean lhsconnected = false;
		Atom lhscand = candidate.getLHS();
		for (AGNode a : agnodes)
		{
			if (lhscand.isNeighbor(a.getLHS()))
			{
				lhsconnected = true;
				break;
			}
		}
		
		boolean rhsconnected = false;
		Atom rhscand = candidate.getRHS();
		for (AGNode a : agnodes)
		{
			if (rhscand.isNeighbor(a.getRHS()))
			{
				rhsconnected = true;
				break;
			}
		}
		
		return (lhsconnected && rhsconnected);
	}
}
