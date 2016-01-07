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

/*
 * A single mapping entity.
 * 
 * Maps a set of atoms to a set of atoms.
 * 
 */
public class HashMapping implements Comparable<HashMapping>, Mapping
{
	
	private HashMap<Atom, Atom> mapping = new HashMap<Atom, Atom>();
	private int gcost;
	private int hcost;
	private int fcost;
	private double featdiff;
	private int hash = 0;
	
	public static int[] PRIMES = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,1601,1607,1609,1613,1619,1621,1627,1637,1657,1663,1667,1669,1693,1697,1699,1709,1721,1723,1733,1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,1823,1831,1847,1861,1867,1871,1873,1877,1879,1889,1901,1907,1913,1931,1933,1949,1951,1973,1979,1987,1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,2063,2069,2081,2083,2087,2089,2099,2111,2113,2129,2131,2137,2141,2143,2153,2161,2179,2203,2207,2213,2221,2237,2239,2243,2251,2267,2269,2273,2281,2287,2293,2297,2309,2311,2333,2339,2341,2347,2351,2357,2371,2377,2381,2383,2389,2393,2399,2411,2417,2423,2437,2441,2447,2459,2467,2473,2477,2503,2521,2531,2539,2543,2549,2551,2557,2579,2591,2593,2609,2617,2621,2633,2647,2657,2659,2663,2671,2677,2683,2687,2689,2693,2699,2707,2711,2713,2719,2729,2731,2741,2749,2753,2767,2777,2789,2791,2797,2801,2803,2819,2833,2837,2843,2851,2857,2861,2879,2887,2897,2903,2909,2917,2927,2939,2953,2957,2963,2969,2971,2999,3001,3011,3019,3023,3037,3041,3049,3061,3067,3079,3083,3089,3109,3119,3121,3137,3163,3167,3169,3181,3187,3191,3203,3209,3217,3221,3229,3251,3253,3257,3259,3271,3299,3301,3307,3313,3319,3323,3329,3331,3343,3347,3359,3361,3371,3373,3389,3391,3407,3413,3433,3449,3457,3461,3463,3467,3469,3491,3499,3511,3517,3527,3529,3533,3539,3541,3547,3557,3559,3571};
	
	// new empty mapping
	public HashMapping()
	{
		
	}

	public HashMapping(HashMap<Atom, Atom> m)
	{
		this.mapping = new HashMap<Atom, Atom>(m);
	}

	public HashMapping clone()
	{
		HashMapping m = new HashMapping();
		m.gcost = this.gcost;
		m.fcost = this.fcost;
		m.hcost = this.hcost;
		m.featdiff = this.featdiff;
		
		for (Atom a : this.mapping.keySet())
		{
			m.extend(a, this.mapping.get(a));
		}
		return m;
	}

	public HashMap<Atom, Atom> getMap()
	{
		return mapping;
	}
	
	public boolean contains(Atom rhs)
	{
		return mapping.containsValue(rhs);
	}
	
	public boolean equals(Object ob)
	{
		HashMapping other = (HashMapping)ob;

		// go over the key-value pairs, and check if the same pair exists in the other HashMap
		Iterator<Atom> it = mapping.keySet().iterator(); 
		while(it.hasNext()) 
		{
			Atom left = (Atom)it.next(); 
			Atom right = (Atom)mapping.get(left);
			
			if (other.getImage(left) != right)
				return false;
		}

		// all this.mapping pairs were in other HashMap, thus if they are of same size, they are equal
		if (this.getDomain().size() == other.getDomain().size())
			return true;
		
		return false;
	}
	
	// check whether can do
	public void extend(Atom lhs, Atom rhs)
	{
//		this.hashCode += "" + a1.getId() + a2.getId();
		this.mapping.put(lhs, rhs);
	}

	public void remove(Atom lhs)
	{
		this.mapping.remove((Object)lhs);
	}
	
	public Collection<Atom> getDomain()
	{
		return this.mapping.keySet();
	}

	public Atom getImage(Atom a)
	{
		return mapping.get(a);
	}


	public Collection<Atom> getRange()
	{
		return this.mapping.values();
	}

	public int hashCode()
	{
		// hashcode never computed, compute it
		if (hash == 0)
		{
			int i = 0;
			long temp = 0;
			
			for (Atom a : mapping.keySet())
			{
				Atom b = mapping.get(a);
				temp += PRIMES[i++] * (a.hashCode()+1) * (b.hashCode()+1);
			}
			
			hash = (int)(temp % Integer.MAX_VALUE);
		}

		return hash;
	}

	public String printMapping()
	{
		StringBuffer sbf = new StringBuffer();
		
		Collection<Atom> ats = this.getDomain();
		List<Atom> ats_list = new ArrayList<Atom>(ats);
		Collections.sort(ats_list);
		sbf.append("ATOM\tROW\tMOL\tATOM\tROW\tMOL\t\n");
		for (Atom a : ats_list)
		{
			sbf.append(a.getSymbol() + "\t" + (a.getMolecularId()+1) + "\t" + (a.getParent().getMolNum()+1));
			sbf.append("\t");
			a = this.getImage(a);
			sbf.append(a.getSymbol() + "\t" + (a.getMolecularId()+1) + "\t" + (a.getParent().getMolNum()+1)
					+ "\n");
		}
		return sbf.toString();
	}

	public void setFCost(int c)
	{
		this.fcost = c;
	}
	
	public void setFCost()
	{
		this.fcost = this.gcost + this.hcost;
	}
	
	public int getFCost()
	{
		return this.fcost;
	}
	
	public void setGCost(int c)
	{
		this.gcost = c;
	}
	
	public int getGCost()
	{
		return this.gcost;
	}
	
	public void setHCost(int c)
	{
		this.hcost = c;
	}
	
	public int getHCost()
	{
		return this.hcost;
	}
	
	public void setCost(int c)
	{
		this.fcost = c;
	}
	
	public int getCost()
	{
		return this.fcost;
	}
	
	public double getFeatDiff()
	{
		return this.featdiff;
	}
	
	public void setFeatDiff(double c)
	{
		this.featdiff = c;
	}
	
	public int size()
	{
		return this.mapping.keySet().size();
	}

	public String toString()
	{
		String ret = "";
		for (Atom a : this.mapping.keySet())
		{
			ret += a.getMolecularId() + " -> " + mapping.get(a).getMolecularId() + "\n";
		}
		return ret;
	}

	public int compareTo(HashMapping x)
	{
		if (this.fcost < x.fcost)
			return -1;
		else if (this.fcost > x.fcost)
			return 1;
		return 0;
	}
	
	public boolean connected()
	{
		return pieces().size() == 1;
	}
	
	public Collection<HashMapping> pieces()
	{
		Collection<HashMapping> pieces = new ArrayList<HashMapping>();
		
		HashMapping rest = this.clone();
		
		while (rest.size() > 0)
		{
			Pair<HashMapping,HashMapping> param = findPiece(rest.getDomain().iterator().next());
			HashMapping connected = param.first;
			rest = param.second;
			
			pieces.add(connected);
		}
		
		return pieces;		
	}
	
	private Pair<HashMapping,HashMapping> findPiece(Atom start)
	{
		HashMapping conn = new HashMapping();
		HashMapping rest = this;
		
		conn.extend(start, this.getImage(start));
		rest.remove(start);
		
base:	while (rest.size() > 0)
		{
 	   		for (Atom cand : rest.getDomain())
			{
				for (Atom in : conn.getDomain())
				{
					if (in.isNeighbor(cand))
					{
						conn.extend( cand, getImage(cand) );
						rest.remove( cand );
						continue base;
					}
				}
			}
 	 		
 	   		break;
		}
		
		return new Pair<HashMapping,HashMapping>(conn,rest);
	}
	
	public int costFunction(Atom lhs, Atom rhs)
	{
		Collection<Atom> msa = this.getDomain();

		int mismatch = 0;
//		double featMismatch = 0.0;
		
		for (Atom a1 : msa)
		{
			Bond sub = a1.getBond(lhs);
			Bond prod = this.getImage(a1).getBond(rhs);

			if (sub == null && prod != null) // bond addition
				mismatch++;
			if (sub != null && prod == null) // bond removal
				mismatch++;
			if (GlobalOptions.bc == true)
				if (sub != null && prod != null && sub.getType() != prod.getType()) // bond type change
					mismatch++;
			// otherwise matches
		}

		return mismatch + this.getGCost();
	}


	// Compute the g(x) (the cost so far)
	//  = the number of differences in the mapped bonds (addition, removal, type change)
	public int costFunction()
	{
		Collection<Atom> msa = this.getDomain();
	
		int mismatch = 0;
	
		for (Atom a1 : msa)
		{
			for (Atom a2 : msa)
			{
				Bond sub = a1.getBond(a2);
				Bond prod = this.getImage(a1).getBond(this.getImage(a2));
	
				if (sub == null && prod != null) // bond addition
					mismatch++;
				if (sub != null && prod == null) // bond removal
					mismatch++;
				if (GlobalOptions.bc == true)
					if (sub != null && prod != null && sub.getType() != prod.getType()) // bond type change
						mismatch++;
				// otherwise matches
			}
		}
	
		return mismatch / 2;
	}

	final class Pair<L,R> {
		public L first;
		public R second;
		public Pair(L first, R second){
			this.first = first;
			this.second = second;
		}
	}
}


