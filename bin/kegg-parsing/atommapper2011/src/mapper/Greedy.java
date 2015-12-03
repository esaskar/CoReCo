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

public class Greedy
{
	public static Mapping greedy(Reaction r)
	{
		return freeGreedy(r, new HashMapping(), Integer.MAX_VALUE, GlobalOptions.plusplus);
	}

	public static Mapping greedy(Reaction r, Mapping map)
	{
		return freeGreedy(r, map, Integer.MAX_VALUE, GlobalOptions.plusplus);
	}
	
	// strict greedy: lhs order is fixed, rhs is argmax
	// time complexity: O(n^2)
	public static Mapping strictGreedy(Reaction r, AstarMapping m, /*Atom[] subatoms, Atom[] prodatoms,*/ PriorityArray Q, double upperbound)
	{
		AstarMapping currentmap = m;
		int firstpos = currentmap.size();
		
		// take each lhs-atom in order and try to map it with all rhs-atoms
		for (int i = firstpos; i < Astar.subatoms.length; i++)
		{
			currentmap = new AstarMapping(currentmap); // create a new currentmap out of previous one
			
			Atom umsa = Astar.subatoms[i];
			
			int k = 1;
			Atom min_atom = null;
			int min_cost = Integer.MAX_VALUE / 2;
			double min_featcost = Integer.MAX_VALUE / 2;
			
			// find cheapest product atom
			for (int j = 0; j < Astar.prodatombins.length; j++)
			{
				if (!(currentmap.getCount(j) < Astar.prodatombins[j].length))
					continue;
				
				Atom umpa = Astar.prodatombins[j][currentmap.getCount(j)];
				
				if (/*currentmap.contains(umpa) ||*/ !umsa.getSymbol().equals(umpa.getSymbol()))
					continue;
				
				// compute new gcost without actually adding 'umpa' to mapping
				int gcost = currentmap.costFunction(umsa, umpa);
				int hcost = currentmap.costEstimate(umsa, umpa);
				
				double featcost = 0.0;
				if (GlobalOptions.plusplus)
					featcost = currentmap.getFeatDiff() + r.normalizedFeatCost(umsa, umpa);

				int fcost = gcost + hcost;
				
				// break out of loop if umpa doesn't increase cost == best solution
				if (fcost == currentmap.getFCost() && featcost <= currentmap.getFeatDiff())
				{
					min_atom = umpa;
					break;
				}
				
				if ((fcost < min_cost) ||
					(fcost == min_cost && featcost < min_featcost) ||
					(fcost == min_cost && Math.random() <= ((1/k) + 0.01)))
				{
					min_cost = fcost;
					min_featcost = featcost;
					min_atom = umpa;
					k++;
				}
			}

			// extend: update bd, update gcost, update hcost
			currentmap.extend(umsa, min_atom);
			
			if (GlobalOptions.plusplus)
				currentmap.setFeatDiff(currentmap.getFeatDiff() + r.normalizedFeatCost(umsa, min_atom));
			
			// return partial mapping if no cheap enough is found
			if (currentmap.getFCost() > upperbound || (GlobalOptions.one && currentmap.getFCost() == upperbound))
				break;
			
			if (i > firstpos) // don't add first greedy mapping, that's already in Q
				Q.offer(currentmap);
		}
		
		return currentmap;
	}
	
	// free greedy algorithm: chooses from the set (lhsborder X rhsborder) the minimum one
	// time complexity: O(n^3)
	public static Mapping freeGreedy(Reaction r, Mapping m, double upperbound, boolean useFeats)
	{
		Collection<Atom> UMPA = r.getProdsAtoms();
		UMPA.removeAll(m.getRange());
		
		Atom[] atoms = new Atom[r.getSubsAtoms().size()];
		for (Atom a : r.getSubsAtoms())
			atoms[a.getOrder()] = a;
		
		
		for (int i = 0; i < atoms.length; i++)
		{
			Atom lhs = atoms[i];
			
			Atom best = null;
			double min_cost = Integer.MAX_VALUE;
			int min_bondcost = Integer.MAX_VALUE;
			
			for (Atom rhs : UMPA)
			{
				if (!lhs.getSymbol().equals(rhs.getSymbol()))
					continue;

				int bondcost = m.costFunction(lhs, rhs);

				if (useFeats)
				{
					double featcost = r.normalizedFeatCost(lhs, rhs);

					if (featcost < min_cost ||
							(featcost <= min_cost && bondcost < min_bondcost))
					{
						best = rhs;
						min_cost = featcost;
						min_bondcost = bondcost;
					}
				}
				else
				{
					if (bondcost < min_bondcost)
					{
						best = rhs;
						min_bondcost = bondcost;
					}
				}
			}

			UMPA.remove(best);

			m.extend(lhs, best);
			m.setGCost(min_bondcost);
//			m.setFCost();
		}
		
		return m;
	}
}


