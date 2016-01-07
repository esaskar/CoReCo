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

public class DFS
{
	private static Collection<Mapping> result;
	private static double ub;

	public static Collection<Mapping> runDFS(Reaction r)
	{
		return runDFS(r, new HashMapping(), Integer.MAX_VALUE, GlobalOptions.plusplus);
	}
	
	public static Collection<Mapping> runDFS(Reaction r, Mapping m, double ub, boolean useFeats)
	{
		result = new HashSet<Mapping>();
		DFS.ub = ub;

		Collection<Atom> UMSA = r.getSubsAtoms();
		UMSA.removeAll(m.getDomain());
		Collection<Atom> UMPA = r.getProdsAtoms();
		UMPA.removeAll(m.getRange());
		
		dfs(r, m, new TreeSet<Atom>(UMSA), new LinkedList<Atom>(UMPA), useFeats);

		double min_cost = Integer.MAX_VALUE;

		// remove too large costs
		for (Mapping ma : result)
		{
			if (ma.getCost() < min_cost)
				min_cost = ma.getCost();
		}

		Iterator<Mapping> it = result.iterator();
		while (it.hasNext())
		{
			if (it.next().getCost() > min_cost)
				it.remove();
		}

		return result;
	}

	private static void dfs(Reaction r, Mapping m, SortedSet<Atom> UMSA,
			LinkedList<Atom> UMPA, boolean useFeats)
	{
		// mapping complete (we are at leaf level)
		if (m.getDomain().size() == r.getSubsAtoms().size())
		{
			int cost = m.costFunction();

			if (cost <= ub)
			{
				ub = cost;
				Mapping m_new = m.clone();
				m_new.setGCost(cost);
				result.add(m_new);
			}

			return;
		}

		// not at leaf level, select alphabetically first atom and try to pair it
		// with all remaining nodes in RHS
		Atom subatom = UMSA.first();
		UMSA.remove(subatom);

		for (int i = 0; i < UMPA.size(); i++)
		{
			Atom umpa = UMPA.removeFirst();

			if (umpa.getSymbol().equals(subatom.getSymbol()))
			{
				double cost = m.costFunction(subatom, umpa);
				
				m.extend(subatom, umpa);

				if (cost <= ub)
					dfs(r, m, UMSA, UMPA, useFeats);

				m.remove(subatom);
			}

			UMPA.addLast(umpa);
		}

		UMSA.add(subatom);
	}
}
