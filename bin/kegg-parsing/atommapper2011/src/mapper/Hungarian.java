
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
/**
 * 
 * Created by: Markus Heinonen
 * 
 */
package mapper;

import java.util.*;
import java.io.*;

public class Hungarian
{
	// mate of both sides of nodes, |V| + |U| (2*n size)
	private static int mate[];
	// label, |V| + |U| (2*n size)
	private static int label[];
	// slack[j] is the minimum alpha for j
	private static int slack[];
	private static final int MAX_VALUE = Integer.MAX_VALUE;
	private static final int MIN_VALUE = Integer.MIN_VALUE;
	private static final int NO_PARENT = -2;
	private static final int NONE = -1;
	// size of graphs, total 2*n
	private static int n;
	private static int matepairs;
	// The nodes of V of alternative tree
	//  indicates the parent in T
	private static int S[];
	private static int Scount;
	// for each node, whether it is in S or not
	private static boolean Sb[];
	// Indicates the parent in S
	private static int T[];
	private static int Tcount;
	private static boolean Tb[];
	// Lists neighbors of S
	private static boolean N[];
	private static int Ncount;
	private static int[][] weights;
	private static int root;

	public static List<Mapping> match(Reaction r)
	{
		Collection<Atom> subs = r.getSubsAtoms();
		Collection<Atom> prods = r.getProdsAtoms();

		int[][] weights = new int[subs.size()][prods.size()];
		
		for (Atom sub_atom : subs)
		{
			for (Atom prod_atom : prods)
			{
				if (GlobalOptions.plusplus)
				{
					// use features to compute the similarity

					int i = sub_atom.getOrder();
					int j = prod_atom.getOrder();
					
					if ( !sub_atom.getSymbol().equals(prod_atom.getSymbol()))
					{
						weights[i][j] = 0;
						continue;
					}
					
					weights[i][j] = (int)(100 * (1.0 - r.normalizedFeatCost(sub_atom, prod_atom)));
				}
				
				else
				{
					// place weights to atom-pairs without feat-information
					//  atom matches -> 30% base score
					//  no atom match -> 0% score
					//  atomneighbors match -> 50% score points
					//  bondneighbors match -> 20% score points
					// 
					// thus the formula is 
					//
					// 30 + (match_atompairs / all_atompairs) * 50
					// + (match_bondpairs / all_bondpairs) * 20
					//
									
					int i = sub_atom.getOrder();
					int j = prod_atom.getOrder();
					
					Molecule subneigh = new Molecule(sub_atom.getAtomNeighbors(), sub_atom.getBondNeighbors(), false);
					Molecule prodneigh = new Molecule(prod_atom.getAtomNeighbors(), prod_atom.getBondNeighbors(), false);
					
					if ( !sub_atom.getSymbol().equals(prod_atom.getSymbol()))
					{
						weights[i][j] = 0;
						continue;
					}
					
					weights[i][j] = 30;
					
					double res1 = subneigh.atomSpectrumSimilarity(prodneigh);
					double res2 = subneigh.bondSpectrumSimilarity(prodneigh);
					assert (res1 >= 0 && res2 >= 0);
					weights[i][j] += (int)(res1 * 50);
					weights[i][j] += (int)(res2 * 20);
				}
			}
		}
		
		int[] mapping = match(weights);
		
		Mapping m = new HashMapping();
		
		for (int i = 0; i < mapping.length/2; i++)
		{
			for (Atom s : subs)
			{
				if (s.getOrder() == i)
				{
					for (Atom p : prods)
					{
						if (p.getOrder() == mapping[i])
						{
							m.extend(s, p);
							break;
						}
					}
				}
			}
		}
		
		m.setGCost(m.costFunction());
//		m.setFCost();
		
		ArrayList<Mapping> ms = new ArrayList<Mapping>();
		ms.add(m);
		System.out.println("Score " + computeScore(weights, mapping) + " out of " + subs.size()*100);
		
		return ms;
	}
	
	public static int[] match(int weights[][])
	{
		n = weights.length;
		Hungarian.weights = weights;

		if (n > 0)
		{
			// |V| + |U| (nodeptr >= 0, none = -1)
			mate = new int[2 * n];
			// |V| + |U| (integer >= 0, none = -1)
			label = new int[2 * n];
			matepairs = 0;
			S = new int[n];
			Sb = new boolean[n];
			Scount = 0;
			T = new int[n];
			Tb = new boolean[n];
			Tcount = 0;
			N = new boolean[n];
			Ncount = 0;
			slack = new int[n];
			root = -1;
		}

		// set mates none
		for (int i = 0; i < 2 * n; i++)
			mate[i] = -1;

		// (1.a) Generate initial labeling
		// - set U-labels zero, set V-labels max of edge-weights

		// set labels of U to zero
		for (int j = 0; j < n; j++)
			label[j + n] = 0;

		// set labels of V to max of edge-weights
		for (int i = 0; i < n; i++)
		{
			int max = MIN_VALUE;

			for (int j = 0; j < n; j++)
			{
				if (weights[i][j] > max)
					max = weights[i][j];
			}

			label[i] = max;
		}

		// (1.b) Generate initial matching $M$ in $E_l$
		// - naively try to map each node to its first available neighbor,
		//   guaranteed to find at least a matching of size one

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				// if edge exists and other end is free and labels match weight
				if (weights[i][j] > 0 && mate[n + j] == -1
						&& label[i] + label[n + j] == weights[i][j])
				{
					setmate(i, j);
					matepairs++;
					break;
				}
			}
		}

		// (2) If matching is perfect (mate is full), stop
		//     Otherwise, pick free node from V, set S = {v}, T = emptyset

		// initialization done
		stage2: while (true)
		{
			assert (matepairs <= n);

			// if done, finish!
			if (matepairs == n)
				return mate;

			// pick first free node in V
			root = -1;
			for (int i = 0; i < n; i++)
			{
				if (mate[i] == -1)
				{
					root = i;
					break;
				}
			}

			if (root == -1)
				System.out.println("Could not find free vertex v");

			// Set S, T and N empty
			for (int j = 0; j < n; j++)
			{
				Tb[j] = false;
				Sb[j] = false;
				S[j] = NONE;
				T[j] = NONE;
				N[j] = false;
			}

			// set S = {v}
			Scount = 1;
			S[root] = NO_PARENT;
			Sb[root] = true;
			Tcount = 0;
			Ncount = 0;
			extendN(root);

			// initialize slacks
			for (int j = 0; j < n; j++)
				slack[j] = label[root] + label[n + j] - weights[root][j];

			// (3) If N_l(S) == T, update labels (thus N_l(S) becomes != T)
			// - Compute alpha_l
			// - Compute new labeling l'

			// test if neighbors of S \subset V (with label[i]+[j] == weight[i][j])
			//   are the same as T

			stage3: while (true)
			{
				boolean n_equals_t = NequalsT();

				// N(S) == T
				if (n_equals_t == true)
				{
					// (3) - find smallest non-zero difference in label and weight

					// Compute alpha as the minimum of slacks of T

					int alpha = MAX_VALUE;

					for (int j = 0; j < n; j++)
					{
						if (N[j] == false && slack[j] < alpha) // node j is in T (N==T)
							alpha = slack[j];
					}

					if (alpha <= 0)
						System.out.println("Alpha is zero, shouldn't probably be");

					// Compute new labeling

					// all nodes in S (subset of V)
					for (int i = 0; i < n; i++)
						if (Sb[i] == true)
							label[i] -= alpha;

					// all nodes in T (subset of U)
					for (int j = 0; j < n; j++)
						if (N[j] == true)
							label[n + j] += alpha;

					// rest nodes are left intact

					// Compute new slacks
					// for all j not in T, slack -= alpha

					for (int j = 0; j < n; j++)
						if (N[j] == false)
							slack[j] -= alpha;

					// make sure that N is updated
					for (int i = 0; i < n; i++)
						if (Sb[i] == true)
							extendN(i);

					// the changes to labels force such that N(S) != T
					n_equals_t = false;
				}

				// (4) If N(S) != T,
				// - pick y from {N(S) - T}
				// - if y is free, augment M with v - y and go to (2)
				// - if y is mated, set S += {mate{y]}, T += {y}, go to (3)

				if (n_equals_t == false)
				{
					// pick first node in U in (N - T)
					for (int j = 0; j < n; j++)
					{
						if (N[j] == true && Tb[j] == false)
						{
							// node is free
							if (mate[n + j] == -1)
							{
								// because node 'j' is free (in U), we have
								// found an alternating path 
								// insert 'j' to T

								// find the parent
								int parent = findParent(j);

								T[j] = parent;
								Tb[j] = true;
								Tcount++;

								augment(j);

								// go to (2)
								continue stage2;
							}
							// node is matched
							// must update all slacks to check if the new node to be added to S is cheaper
							else
							{
								// find the parent
								int parent = findParent(j);

								// set tree structure
								T[j] = parent;
								Tb[j] = true;
								Tcount++;

								S[mate[n + j]] = j;
								Sb[mate[n + j]] = true;
								Scount++;

								extendN(mate[n + j]);

								// update slacks
								// we have a new node in S, namely mate[n+j], try if its cheaper
								int new_in_S = mate[n + j];
								for (int k = 0; k < n; k++)
								{
									int diff = label[new_in_S] + label[n + k]
											- weights[new_in_S][k];
									if (diff < slack[k])
										slack[k] = diff;
								}

								// go to (3)
								continue stage3;
							}
						}
					}
				}
			}
		}
	}

	public static int computeScore(int[][] weights, int[] result)
	{
		int sum = 0;
		
		for (int i = 0; i < result.length / 2; i++)
		{
			sum += weights[i][result[i]];
		}
		
		return sum;
	}
	
	private static int findParent(int j)
	{
		// find a parent to node 'j'

		// If only root in the tree, return it
		if (Scount == 1)
			return root;

		for (int i = 0; i < n; i++)
		{
			if (Sb[i] == true && label[i] + label[n + j] == weights[i][j])
				return i;
		}

		System.out.println("findParent: could not find parent!");
		return -1;
	}

	private static void extendN(int i)
	{
		// As N can only grow in one phase (with expansion of S), we can extend N
		// with each node inserted into S
		//
		// extend N with neighbors of 'i'

		for (int j = 0; j < n; j++)
		{
			// not in N already, and labels match weight
			if (N[j] == false && label[i] + label[j + n] == weights[i][j])
			{
				N[j] = true;
				Ncount++;
			}
		}
	}

	private static void augment(int j)
	{
		// map j and T[j]
		// call augment2( S[T[j]] )
		// map S[T[j]] and T[S[T[j]]]

		setmate(T[j], j);

		if (T[j] == root)
		{
			matepairs++;
			return;
		} else
			augment(S[T[j]]);
	}

	private static boolean NequalsT()
	{
		// N and T has to match in size
		if (Ncount != Tcount)
			return false;

		// Also both have to have same nodes

		for (int i = 0; i < n; i++)
			if (N[i] != Tb[i])
				return false;

		return true;
	}

	private static void setmate(int i, int j)
	{
		mate[i] = j;
		mate[n + j] = i;
	}

	public static int[][] ParseReactionWeights(String reac)
	{
		int weights[][] = null;

		// read in a molecule
		String filepath = "/home/fs/mqheinon/duuni/atomisota/koodi/" + reac + ".w";
		File f = new File(filepath);
		String temp = "";

		Scanner sc = null;
		try
		{
			sc = new Scanner(f);
		} catch (Exception e)
		{
		}

		temp = sc.nextLine();

		int size = Integer.parseInt(temp);

		weights = new int[size][size];

		int i = 0;
		int j = 0;

		while (sc.hasNextLine())
		{
			j = 0;
			temp = sc.nextLine();

			String w[] = temp.split("\t");

			for (String s : w)
			{
				weights[i][j] = Integer.parseInt(s);
				j++;
			}

			i++;
		}

		return weights;
	}
}


