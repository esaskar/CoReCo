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

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.util.*;

public class Astar
{
	//	 function A*(start,goal)
	//     closedset := the empty set                 % The set of nodes already evaluated.
	//     openset := set containing the initial node % The set of tentative nodes to be evaluated.
	//     g_score[start] := 0                        % Distance from start along optimal path.
	//     while openset is not empty
	//         x := the node in openset having the lowest f_score[] value
	//         if x = goal
	//             return reconstruct_path(came_from,goal)
	//         remove x from openset
	//         add x to closedset
	//         foreach y in neighbor_nodes(x)
	//             if y in closedset
	//                 continue
	//             tentative_g_score := g_score[x] + dist_between(x,y)
	//             tentative_is_better := false
	//             if y not in openset
	//                 add y to openset
	//                 h_score[y] := heuristic_estimate_of_distance_to_goal_from(y)
	//                 tentative_is_better := true
	//             elseif tentative_g_score < g_score[y]
	//                 tentative_is_better := true
	//             if tentative_is_better = true
	//                 came_from[y] := x
	//                 g_score[y] := tentative_g_score
	//                 f_score[y] := g_score[y] + h_score[y] % Estimated total distance from start to goal through y.
	//     return failure	

	public static Atom[] subatoms; // ordered based on 'order'
	public static Atom[] prodatoms; // ordered based on 'order'
	public static Atom[][] prodatombins; // prodatoms are grouped into bins (inner arrays) of isomorphic atoms
	public static int reacsize;
	public static int atomstrsize; // number of unique atomstrs (eg. "C|CCO")
	public static int bondstrsize; // number of unique bondstr (eg. "C-C")
	public static Reaction r;

	public static Collection<Mapping> astarFixedOrder(Reaction r)
	{
		long starttime = System.currentTimeMillis();
		long time = System.currentTimeMillis();
		long greedytime = 0;
		long temptime = 0;
		long statustime = time;
		long processednodes = 0;

		// Results
		Collection<Mapping> result = new HashSet<Mapping>();

		// Pointer to the previously processed node
		AstarMapping prevmap;
		// Pointer to current node
		AstarMapping currentmap;
		// Pointer to newly created node
		AstarMapping newmap;

		// Pointer to used mappings (nodes)
//		HashSet<Mapping> used = new HashSet<Mapping>();


		// reaction size: atom count
		Astar.reacsize = r.getSubsAtoms().size();

		// atom-neighbor-strings, e.g. "C|COO"
		r.computeAtomStringIndex();
		r.computeBondStringIndex();

		// reactions' atoms are ordered [0...n],
		// place them in order to subatoms/prodatoms
		Astar.subatoms = new Atom[Astar.reacsize];
		for (Atom temp : r.getSubsAtoms())
			Astar.subatoms[temp.getOrder()] = temp;

		
		Astar.prodatoms = new Atom[Astar.reacsize];
		for (Atom temp : r.getProdsAtoms())
			Astar.prodatoms[temp.getOrder()] = temp;
		
		Astar.prodatombins = r.getProdAtomBins();

		// initialize Q with empty mapping 'zeromap'
		AstarMapping zeromap = new AstarMapping(r);
		//		Mapping greedymap = Greedy.strictGreedy(r, firstnode, subatoms, prodatoms, Q, costEstimates, Integer.MAX_VALUE);

		int upperbound = r.computeInitialUB();
		//		int upperbound = greedymap.costFunction();
		//		int upperbound = Math.max(greedymap.costFunction(), Hungarian.match(r).iterator().next().getFCost());
		int lowerbound = zeromap.getHCost();

//		// PriorityQueue of nodes ( mapping=(short[],bitset,spectra) )
//		PriorityQueue<AstarMapping> Q = new PriorityQueue<AstarMapping>(100,
//				new AstarMappingComparator()); // unprocessed nodes

		PriorityArray Q = new PriorityArray(upperbound);
		
		// Progress of ub/lb
		// Set all points to -1, and initialize with greedy/bpm ub and h(x) lb
		// thus for (lb,ub) = (2,6) we get
		//  lbhistory = [-1, -1,  0, -1, -1, -1]
		//  ubhistory = [-1, -1, -1, -1, -1,  0]

//		int lbhistory[] = new int[upperbound + 1];
//		int ubhistory[] = new int[upperbound + 1];

//		for (int i = 0; i < upperbound + 1; i++)
//		{
//			lbhistory[i] = -1;
//			ubhistory[i] = -1;
//		}
//
//		lbhistory[lowerbound] = 0;
//		ubhistory[upperbound] = 0;

		Q.offer(zeromap);
//		SearchTree tree = new SearchTree(subatoms);
		
		prevmap = zeromap;
		currentmap = zeromap;

		// start the algorithm
		while (!Q.isEmpty())
		{
			currentmap = Q.poll(); // cheapest
			
			processednodes++;

			if (currentmap.getFCost() > upperbound) // quit fast
				continue;

			// let's discard the result if its not strictly better than previous
			if (GlobalOptions.one && currentmap.getFCost() == upperbound)
				continue;

			if (currentmap.size() == Astar.reacsize) // complete mapping?
			{
				if (currentmap.getFCost() < upperbound) // better result
				{
					result.clear();

					upperbound = currentmap.getFCost();

//					if (ubhistory[upperbound] == -1)
//						ubhistory[upperbound] = processednodes;

//					Collection<AstarMapping> removables = new HashSet<AstarMapping>();
//					for (AstarMapping x : Q)
//						if (x.getFCost() > upperbound || (GlobalOptions.one && x.getFCost() == upperbound))
//							removables.add(x);
//					Q.removeAll(removables);
					
					// PriorityArray version
					int ub_to_remove = GlobalOptions.one ? upperbound : upperbound + 1;
					Q.removeByScore(ub_to_remove);
				}

				if (result.size() < GlobalOptions.maxresults)
					result.add(currentmap.clone());
				else if (upperbound == lowerbound || GlobalOptions.lb >= currentmap.getFCost()) // quit if all 100 results are optimal, or we hit lb threshold
					break;

				continue;
			}

			// lowerbound only grows
			if (currentmap.getFCost() > lowerbound)
				lowerbound = currentmap.getFCost();

//			if (lbhistory[lowerbound] == -1)
//				lbhistory[lowerbound] = processednodes;

			if (upperbound == lowerbound && result.size() >= GlobalOptions.maxresults)
				break;

			// take next atom as lhs
			Atom lhs = Astar.subatoms[currentmap.size()];

			// go through prodatoms, discard if already mapped
			for (int i = 0; i < prodatombins.length; i++)
			{
				if (!(currentmap.getCount(i) < Astar.prodatombins[i].length))
					continue;
				
				Atom rhs = Astar.prodatombins[i][currentmap.getCount(i)];

				if (/*currentmap.contains(rhs) ||*/ !lhs.getSymbol().equals(rhs.getSymbol()))
					continue;

				// create new node
				newmap = new AstarMapping(currentmap); // create the new one out of current node
				newmap.extend(lhs, rhs); // extend (updates costs and bd)

				// Don't go there if cost too large
				if (newmap.getFCost() > upperbound)
					continue;
				
				if (GlobalOptions.one && newmap.getFCost() == upperbound)
					continue;

				if (GlobalOptions.plusplus)
					newmap.setFeatDiff(currentmap.getFeatDiff() + r.normalizedFeatCost(lhs, rhs));
//				else
//					newmap.setFeatDiff(0.0);

				Q.offer(newmap);
//				tree.addNode(currentnode, newnode);
			}

			// Greedy search whenever we backtrack in the search tree
			// Approximate test: if currentnode is not 1 larger than prevnode, do the greedy
			if (upperbound > lowerbound && currentmap.size() < prevmap.size()
					&& processednodes % GlobalOptions.greedyFreq == 0)
			{
				temptime = System.currentTimeMillis();
				/*Mapping greedymap =*/ Greedy.strictGreedy(r, prevmap, /*Astar.subatoms, Astar.prodatoms,*/ Q, upperbound);
				greedytime += System.currentTimeMillis() - temptime;
			}
			
			// check quitting conditions every 100ms
			if ((System.currentTimeMillis() - time) > 100)
			{
				time = System.currentTimeMillis();
				
				// runtime condition
				if ((time - starttime) / 1000 > GlobalOptions.maxtime)
				{
					System.out.println("runtime " + GlobalOptions.maxtime + ", quitting..");
					return null;
				}
				
				// heapsize condition
				if (Q.size() > GlobalOptions.MAX_HEAP) // more than 1 million, -> quit
				{
					System.out.println("MAX HEAP " + GlobalOptions.MAX_HEAP + " hit, quitting...");
					return null;
				}
				
				// out of memory
				long usableFreeMemory = Runtime.getRuntime().maxMemory() - Runtime.getRuntime().totalMemory() + Runtime.getRuntime().freeMemory();
				if (usableFreeMemory < 1048576 * 50)
				{
					System.out.println("Less than 50 megs free memory, quitting...");
					return null;
				}
			}
			
			if ((System.currentTimeMillis() - statustime) > 5000)
			{
//				System.gc();
				
				// check for lb-misses and sort everything
//				Q.refresh();
				
				statustime = System.currentTimeMillis();
				
				// compute distribution
//				int[] dist = new int[upperbound+1];
//				for (AstarMapping m : Q)
//					dist[m.getFCost()]++;
//				String d = "dist:";
//				for (int i = lowerbound; i <= upperbound;i++)
//					d += dist[i] + ",";
				
//				MemoryUsage mem = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
				
				long usableFreeMemory = Runtime.getRuntime().maxMemory() - Runtime.getRuntime().totalMemory() + Runtime.getRuntime().freeMemory();
				// print status information
				time = System.currentTimeMillis();
				System.out.println("Qsize: " + Q.size() + "\tRes: " + result.size()
						+ "\tUsed: " + processednodes + "\tlb " + lowerbound + "\tub "
						+ upperbound + "\t@ " + (time - starttime) / 1000 + "s ("
						+ greedytime / 1000 + "s greedy). (free: "
						+ usableFreeMemory / 1048576 + "M)");
//						+ Runtime.getRuntime().freeMemory() / 1048576 + "M, jvm: "
//						+ Runtime.getRuntime().totalMemory() / 1048576 + "M, max: "
//						+ Runtime.getRuntime().maxMemory() / 1048576 + "M)" 
//						+ /*d +*/ " (" + mem.getUsed()/1048576 + "M," + mem.getCommitted()/1048576 + "M," + mem.getMax()/1048576 + "M)");
				System.out.println(Q);


//				prevprocessed = processednodes;
			}

			prevmap = currentmap;

//			tree.drawTree();
		}

		////////////////////////////////////

//		if (lbhistory[lowerbound] == -1)
//			lbhistory[lowerbound] = processednodes;
//		if (ubhistory[upperbound] == -1)
//			ubhistory[upperbound] = processednodes;
//
//		if (GlobalOptions.progressinfo)
//			r.writeProgressInfo(lbhistory, ubhistory);

		return result;
	}
}
