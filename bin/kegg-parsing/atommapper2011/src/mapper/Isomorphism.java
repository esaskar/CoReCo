
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
/*
 * Isomorphism algorithm implementation (VF2) for reaction graphs
 * 
 */
package mapper;

import java.util.*;


public class Isomorphism
{
	class AtomPair
	{
		public RGAtom lhs;
		public RGAtom rhs;
		public AtomPair(RGAtom lhs, RGAtom rhs)
		{
			this.lhs = lhs;
			this.rhs = rhs;
		}
	}
	
	class StateNode
	{
		public Map<RGAtom,RGAtom> map;
		public BitSet lhsborder;
		public BitSet rhsborder;
		
		StateNode()
		{
//			map = new HashMap<RGAtom,RGAtom>();
//			lhsborder = new BitSet();
//			rhsborder = new BitSet();
		}
		
		public int size()
		{
			return map.size();
		}
		
		public String toString()
		{
			String s = "";
			for (RGAtom a : map.keySet())
				s += a + " <-> " + map.get(a) + "\n";
			return s;
		}
	}
	
	


	private ReactionGraph rg1;
	private ReactionGraph rg2;
	private Direction dir;
	private int itercount = 0;
	private boolean iso = false;
	private Collection<Map<RGAtom,RGAtom>> solutions;
	private RGAtom[] g1atoms;
	private RGAtom[] g2atoms;
	private int size;
	
	public Isomorphism(ReactionGraph rg1, ReactionGraph rg2)
	{
		this.rg1 = rg1;
		this.rg2 = rg2;
		solutions = new ArrayList<Map<RGAtom,RGAtom>>();
		iso = false;
		
		g1atoms = new RGAtom[rg1.size()];
		g2atoms = new RGAtom[rg2.size()];
		
		for (RGAtom a : rg1.getAtoms())
			g1atoms[a.getId()] = a;
		for (RGAtom a : rg2.getAtoms())
			g2atoms[a.getId()] = a;
	}
	
	
	public boolean VF2()
	{
		if (this.rg1.size() != this.rg2.size())
			return false;
		
		size = rg1.size();
		
		// otherwise continue with algorithm
		StateNode start = new StateNode();
		start.map = new HashMap<RGAtom,RGAtom>();
		start.lhsborder = new BitSet(size);
		start.rhsborder = new BitSet(size);
		dir = Direction.FORWARD;
		match(start);
		
		if (iso)
			return true;
		
		start = new StateNode();
		start.map = new HashMap<RGAtom,RGAtom>();
		start.lhsborder = new BitSet(size);
		start.rhsborder = new BitSet(size);
		dir = Direction.BACKWARD;
		match(start);
		
		return iso;
	}
	
	private void match(StateNode s)
	{
		if (iso)
			return;
		
		itercount++;
		
		// full isomorphic mapping found
		if (s.size() == size)
		{
			iso = true;
//			solutions.add( new HashMap<RGAtom,RGAtom>(s.map) ); // shallow copy
			return;
		}
		
		for (AtomPair p : candidatepairs(s))
		{
			RGAtom lhs = p.lhs;
			RGAtom rhs = p.rhs;
			
			if (feasible(s, lhs, rhs))
			{
				StateNode ns = new StateNode();
				ns.map = new HashMap<RGAtom,RGAtom>(s.map); // shallow copy, ok because RGAtom's don't ever change
				ns.map.put(lhs, rhs);
				
				ns.lhsborder = (BitSet)s.lhsborder.clone();
				ns.lhsborder.set(lhs.getId(), false);
				ns.rhsborder = (BitSet)s.rhsborder.clone();
				ns.rhsborder.set(rhs.getId(), false);
				
				
//				ns.lhsborder = new ArrayList(s.lhsborder);
////				Collections.copy(ns.lhsborder, s.lhsborder); // shallow copy
//				ns.lhsborder.remove(lhs); // might not contain
//				ns.rhsborder = new ArrayList(s.rhsborder);
////				Collections.copy(ns.rhsborder, s.rhsborder); // shallow copy
//				ns.rhsborder.remove(rhs); // might not contain
				
				for (RGAtom ne : lhs.getAtomNeighbors())
					if (ns.lhsborder.get(ne.getId()) == false && !ns.map.containsKey(ne))
						ns.lhsborder.set(ne.getId());
				for (RGAtom ne : rhs.getAtomNeighbors())
					if (ns.rhsborder.get(ne.getId()) == false && !ns.map.containsValue(ne))
						ns.rhsborder.set(ne.getId());

//				Collections.sort(ns.lhsborder);
//				Collections.sort(ns.rhsborder);
				
				match(ns);
			}
		}
	}
	
	
	private List<AtomPair> candidatepairs(StateNode s)
	{
		List<AtomPair> pairs = new ArrayList<AtomPair>();
		
		// both sides have border to go through
		if (!s.lhsborder.isEmpty() && !s.rhsborder.isEmpty())
		{
			RGAtom rhs = null;
			for (int i = 0; i < s.rhsborder.size(); i++)
			{
				if (s.rhsborder.get(i) == true)
				{
					rhs = g2atoms[i];
					break;
				}
			}
			
//			RGAtom rhs = s.rhsborder.get(0);
			for (int i = 0; i < s.lhsborder.size(); i++)
			{
				if (s.lhsborder.get(i) == true)
					pairs.add( new AtomPair(g1atoms[i],rhs) );
			}
//			for (RGAtom lhs : s.lhsborder)
//			{
//				pairs.add(new AtomPair(lhs,rhs));
//			}
		}
		// either side has emptied its border
		else
		{
			// set of all atoms not mapped on 'rg2'
			Set<RGAtom> rhsleft = new HashSet<RGAtom>(rg2.getAtoms());
			rhsleft.removeAll(s.map.values());
			// pick smallest id atom from them
			int min_id = 100000;
			RGAtom min_atom = null;
			for (RGAtom a : rhsleft)
			{
				if (a.getId() < min_id)
				{
					min_id = a.getId();
					min_atom = a;
				}
			}
		
			for (RGAtom lhs : rg1.getAtoms())
				if (!s.map.containsKey(lhs))
					pairs.add(new AtomPair(lhs,min_atom));
		}
		
		return pairs;
	}
	
	private boolean feasible(StateNode s, RGAtom lhs, RGAtom rhs)
	{
		// the regions of G1 and G2 are divided into three areas:
		// - (MR) mapped region
		// - (BR) adjacent borders to mapped regions
		// - (RR) remote regions
		//
		// the regions don't overlap
		//
		
		// atom symbol has to match
		if (!lhs.getSymbol().equals(rhs.getSymbol()))
			return false;
		
		// check that lhs's mapped neighbors match rhs's mapped neighbors
		Set<RGAtom> lhs_mapped_neighs = new HashSet<RGAtom>();
		for (RGAtom ne : lhs.getAtomNeighbors())
			if (s.map.containsKey(ne))
				lhs_mapped_neighs.add( ne );
		Set<RGAtom> rhs_mapped_neighs = new HashSet<RGAtom>();
		for (RGAtom ne : rhs.getAtomNeighbors())
			if (s.map.containsValue(ne))
				rhs_mapped_neighs.add( ne );
		
		// mapped neighborhoods have to match
		if (lhs_mapped_neighs.size() != rhs_mapped_neighs.size())
			return false;
		
		// check that each atom matches some atom from other side through mapping
		for (RGAtom a : lhs_mapped_neighs)
		{
			if (!rhs_mapped_neighs.contains(s.map.get(a)))
				return false;
			
			if (dir == Direction.FORWARD       && a.getBond(lhs).getChangetype() != s.map.get(a).getBond(rhs).getChangetype() )
				return false;
			else if (dir == Direction.BACKWARD && a.getBond(lhs).getChangetype() != -1*s.map.get(a).getBond(rhs).getChangetype() )
				return false;
		}
		
		
		// border area sizes have to match
		int lhsneighs = 0;
		for (RGAtom ne : lhs.getAtomNeighbors())
			if (s.lhsborder.get(ne.getId()) == true)
				lhsneighs++;
		int rhsneighs = 0;
		for (RGAtom ne : rhs.getAtomNeighbors())
			if (s.rhsborder.get(ne.getId()) == true)
				rhsneighs++;
		
		if (lhsneighs != rhsneighs)
			return false;
		
		// remote area sizes have to match
		int lhsremotes = 0;
		for (RGAtom ne : lhs.getAtomNeighbors())
			if (s.lhsborder.get(ne.getId()) == false && !s.map.containsKey(ne))
				lhsremotes++;
		int rhsremotes = 0;
		for (RGAtom ne : rhs.getAtomNeighbors())
			if (s.rhsborder.get(ne.getId()) == false && !s.map.containsValue(ne))
				rhsremotes++;
		
		if (lhsremotes != rhsremotes)
			return false;
		
		return true;
	}
}
